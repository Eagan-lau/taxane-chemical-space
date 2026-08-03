from __future__ import annotations

import csv
import json
import sqlite3
import time
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .chemistry import (
    atom_delta_matches,
    formula_delta,
    parse_reaction_delta_fingerprint,
    require_rdkit,
    stable_hash,
)
from .io import ensure_dir, read_table, write_json
from .screen import _bit_vector_as_int


ALLOWED_ATOMIC_NUMBERS = {
    1,  # H
    5,  # B
    6,  # C
    7,  # N
    8,  # O
    9,  # F
    14,  # Si
    15,  # P
    16,  # S
    17,  # Cl
    34,  # Se
    35,  # Br
    53,  # I
}


@dataclass(frozen=True)
class CompiledGrammarRule:
    grammar_rule_id: str
    smarts_rule_id: str
    semantic_group_id: str
    reaction_type: str
    reaction_smarts: str
    evidence_layer: str
    final_rule_confidence: str
    expected_delta_text: str
    expected_delta: dict[str, int]
    query: Any
    reaction: Any
    query_fp: int


def _text(record: dict[str, Any], key: str) -> str:
    value = record.get(key, "")
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value


def _compile_rules(grammar: pd.DataFrame) -> tuple[list[CompiledGrammarRule], list[dict[str, str]]]:
    Chem, AllChem, *_rest = require_rdkit()
    compiled: list[CompiledGrammarRule] = []
    failures: list[dict[str, str]] = []
    for record in grammar.fillna("").to_dict("records"):
        smarts = _text(record, "reaction_smarts")
        rule_id = _text(record, "grammar_rule_id") or _text(record, "smarts_rule_id")
        try:
            reaction = AllChem.ReactionFromSmarts(smarts)
            if reaction is None:
                raise ValueError("ReactionFromSmarts returned None")
            reaction.Initialize()
            if reaction.GetNumReactantTemplates() != 1:
                raise ValueError("generative grammar requires one reactant template")
            if reaction.GetNumProductTemplates() != 1:
                raise ValueError("generative grammar requires one product template")
            query = reaction.GetReactantTemplate(0)
            query_fp = _bit_vector_as_int(Chem.PatternFingerprint(query, fpSize=2048))
            expected_text = (
                _text(record, "structural_element_delta")
                or _text(record, "effective_reaction_delta")
                or _text(record, "reaction_delta_fingerprint")
            )
            expected = parse_reaction_delta_fingerprint(expected_text)
            compiled.append(
                CompiledGrammarRule(
                    grammar_rule_id=rule_id,
                    smarts_rule_id=_text(record, "smarts_rule_id"),
                    semantic_group_id=_text(record, "semantic_group_id"),
                    reaction_type=_text(record, "reaction_type"),
                    reaction_smarts=smarts,
                    evidence_layer=(
                        _text(record, "evidence_layer_best")
                        or _text(record, "exclusive_release_tier")
                    ),
                    final_rule_confidence=_text(record, "final_rule_confidence"),
                    expected_delta_text=expected_text,
                    expected_delta=expected,
                    query=query,
                    reaction=reaction,
                    query_fp=query_fp,
                )
            )
        except Exception as exc:
            failures.append(
                {
                    "grammar_rule_id": rule_id,
                    "smarts_rule_id": _text(record, "smarts_rule_id"),
                    "reaction_smarts": smarts,
                    "compile_error": str(exc),
                }
            )
    return compiled, failures


def _clear_atom_maps(mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def _origin_atom_indices(product) -> set[int]:
    result: set[int] = set()
    for atom in product.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        if atom.HasProp("react_atom_idx"):
            try:
                result.add(int(atom.GetProp("react_atom_idx")))
            except (TypeError, ValueError):
                pass
    return result


def _origin_retention(product, source) -> float:
    source_heavy = {
        atom.GetIdx() for atom in source.GetAtoms() if atom.GetAtomicNum() > 1
    }
    if not source_heavy:
        return 0.0
    return len(_origin_atom_indices(product) & source_heavy) / len(source_heavy)


def _observed_changed_source_atoms(source, product) -> int:
    product_by_origin: dict[int, Any] = {}
    for atom in product.GetAtoms():
        if atom.HasProp("react_atom_idx"):
            try:
                product_by_origin[int(atom.GetProp("react_atom_idx"))] = atom
            except (TypeError, ValueError):
                pass
    changed = 0
    for source_atom in source.GetAtoms():
        product_atom = product_by_origin.get(source_atom.GetIdx())
        if product_atom is None:
            changed += 1
            continue
        source_state = (
            source_atom.GetAtomicNum(),
            source_atom.GetFormalCharge(),
            source_atom.GetDegree(),
            source_atom.GetIsAromatic(),
        )
        product_state = (
            product_atom.GetAtomicNum(),
            product_atom.GetFormalCharge(),
            product_atom.GetDegree(),
            product_atom.GetIsAromatic(),
        )
        if source_state != product_state:
            changed += 1
            continue
        source_origin_neighbors = sorted(
            neighbor.GetIdx() for neighbor in source_atom.GetNeighbors()
        )
        product_origin_neighbors = sorted(
            int(neighbor.GetProp("react_atom_idx"))
            for neighbor in product_atom.GetNeighbors()
            if neighbor.HasProp("react_atom_idx")
        )
        if source_origin_neighbors != product_origin_neighbors:
            changed += 1
    return changed


def _molecule_identity(mol) -> dict[str, Any]:
    Chem, _AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors, _standard = (
        require_rdkit()
    )
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    full_key = Chem.MolToInchiKey(mol)
    if not full_key:
        raise ValueError("empty InChIKey")
    return {
        "smiles": canonical,
        "full_inchikey": full_key,
        "connectivity_key": full_key[:14],
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": round(rdMolDescriptors.CalcExactMolWt(mol), 6),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "clogp": round(Crippen.MolLogP(mol), 6),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 6),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 6),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
    }


@lru_cache(maxsize=1)
def _morgan_generator():
    from rdkit.Chem import rdFingerprintGenerator

    return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _fingerprint(mol):
    return _morgan_generator().GetFingerprint(mol)


def _initialize_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path, timeout=120.0)
    connection.execute("PRAGMA busy_timeout=120000")
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS nodes (
            space_id TEXT PRIMARY KEY,
            generation_first INTEGER NOT NULL,
            source_node_ids TEXT NOT NULL,
            molecule_names TEXT NOT NULL,
            smiles TEXT NOT NULL,
            full_inchikey TEXT NOT NULL UNIQUE,
            connectivity_key TEXT NOT NULL,
            formula TEXT NOT NULL,
            exact_mass REAL NOT NULL,
            heavy_atom_count INTEGER NOT NULL,
            clogp REAL NOT NULL,
            tpsa REAL NOT NULL,
            hbd INTEGER NOT NULL,
            hba INTEGER NOT NULL,
            rotatable_bonds INTEGER NOT NULL,
            ring_count INTEGER NOT NULL,
            fraction_csp3 REAL NOT NULL,
            formal_charge INTEGER NOT NULL,
            known_g0_full_match INTEGER NOT NULL,
            known_g0_connectivity_match INTEGER NOT NULL,
            known_g0_match_ids TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_nodes_generation
          ON nodes(generation_first);
        CREATE INDEX IF NOT EXISTS idx_nodes_connectivity
          ON nodes(connectivity_key);

        CREATE TABLE IF NOT EXISTS derivation_events (
            event_id INTEGER PRIMARY KEY AUTOINCREMENT,
            generation INTEGER NOT NULL,
            source_space_id TEXT NOT NULL,
            target_space_id TEXT NOT NULL,
            target_is_new INTEGER NOT NULL,
            target_generation_first INTEGER NOT NULL,
            grammar_rule_id TEXT NOT NULL,
            smarts_rule_id TEXT NOT NULL,
            semantic_group_id TEXT NOT NULL,
            reaction_type TEXT NOT NULL,
            evidence_layer TEXT NOT NULL,
            final_rule_confidence TEXT NOT NULL,
            expected_element_delta TEXT NOT NULL,
            observed_element_delta TEXT NOT NULL,
            source_atom_retention REAL NOT NULL,
            observed_changed_source_atoms INTEGER NOT NULL,
            source_product_tanimoto REAL NOT NULL,
            raw_product_index INTEGER NOT NULL,
            known_g0_full_match INTEGER NOT NULL,
            known_g0_connectivity_match INTEGER NOT NULL,
            immediate_reverse_cycle INTEGER NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_events_generation
          ON derivation_events(generation);
        CREATE INDEX IF NOT EXISTS idx_events_source
          ON derivation_events(source_space_id);
        CREATE INDEX IF NOT EXISTS idx_events_target
          ON derivation_events(target_space_id);

        CREATE TABLE IF NOT EXISTS application_audit (
            generation INTEGER NOT NULL,
            source_space_id TEXT NOT NULL,
            grammar_rule_id TEXT NOT NULL,
            source_site_match_count INTEGER NOT NULL,
            raw_product_tuple_count INTEGER NOT NULL,
            accepted_event_count INTEGER NOT NULL,
            unique_accepted_product_count INTEGER NOT NULL,
            rejected_product_count INTEGER NOT NULL,
            application_status TEXT NOT NULL,
            rejection_reason_counts TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS rejection_events (
            generation INTEGER NOT NULL,
            source_space_id TEXT NOT NULL,
            grammar_rule_id TEXT NOT NULL,
            raw_product_index INTEGER NOT NULL,
            rejection_reason TEXT NOT NULL,
            product_smiles TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS generation_parent_progress (
            generation INTEGER NOT NULL,
            source_space_id TEXT NOT NULL,
            parent_ordinal INTEGER NOT NULL,
            status TEXT NOT NULL,
            PRIMARY KEY (generation, source_space_id)
        );
        CREATE INDEX IF NOT EXISTS idx_parent_progress_generation_ordinal
          ON generation_parent_progress(generation, parent_ordinal);
        """
    )
    connection.commit()
    return connection


def _export_query(connection: sqlite3.Connection, query: str, path: Path) -> None:
    cursor = connection.execute(query)
    columns = [description[0] for description in cursor.description]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        while True:
            rows = cursor.fetchmany(10_000)
            if not rows:
                break
            writer.writerows(rows)


def _load_g0(
    connection: sqlite3.Connection,
    nodes_path: Path,
) -> tuple[dict[str, str], dict[str, set[str]], dict[str, Any]]:
    Chem, *_rest = require_rdkit()
    raw = read_table(nodes_path).fillna("")
    smiles_column = next(
        column
        for column in ("standardized_smiles", "smiles", "smiles_raw")
        if column in raw.columns
    )
    id_column = next(
        column
        for column in ("molecule_id", "node_id", "source_id", "row_index")
        if column in raw.columns
    )
    name_column = next(
        (
            column
            for column in ("molecule_name", "input_molecule_name", "name")
            if column in raw.columns
        ),
        None,
    )
    grouped: dict[str, dict[str, Any]] = {}
    invalid = 0
    for record in raw.to_dict("records"):
        mol = Chem.MolFromSmiles(_text(record, smiles_column))
        if mol is None:
            invalid += 1
            continue
        identity = _molecule_identity(mol)
        group = grouped.setdefault(
            identity["full_inchikey"],
            {
                "identity": identity,
                "source_ids": [],
                "names": [],
            },
        )
        group["source_ids"].append(_text(record, id_column))
        if name_column and _text(record, name_column):
            group["names"].append(_text(record, name_column))

    full_to_space: dict[str, str] = {}
    connectivity_to_ids: dict[str, set[str]] = {}
    molecules: dict[str, Any] = {}
    for index, (full_key, group) in enumerate(sorted(grouped.items()), start=1):
        identity = group["identity"]
        space_id = f"G0_{index:05d}"
        source_ids = ";".join(sorted(set(group["source_ids"])))
        names = ";".join(sorted(set(group["names"])))
        connection.execute(
            """
            INSERT INTO nodes VALUES (
              ?, 0, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, 1, ?
            )
            """,
            (
                space_id,
                source_ids,
                names,
                identity["smiles"],
                full_key,
                identity["connectivity_key"],
                identity["formula"],
                identity["exact_mass"],
                identity["heavy_atom_count"],
                identity["clogp"],
                identity["tpsa"],
                identity["hbd"],
                identity["hba"],
                identity["rotatable_bonds"],
                identity["ring_count"],
                identity["fraction_csp3"],
                identity["formal_charge"],
                source_ids,
            ),
        )
        full_to_space[full_key] = space_id
        connectivity_to_ids.setdefault(identity["connectivity_key"], set()).update(
            group["source_ids"]
        )
        molecules[space_id] = mol
    connection.commit()
    return full_to_space, connectivity_to_ids, {
        "input_records": int(len(raw)),
        "valid_records": int(sum(len(group["source_ids"]) for group in grouped.values())),
        "invalid_records": invalid,
        "unique_full_stereo_structures": int(len(grouped)),
    }


def _generation_parents(
    connection: sqlite3.Connection,
    generation: int,
) -> Iterable[tuple[str, str, str]]:
    yield from connection.execute(
        """
        SELECT space_id, smiles, full_inchikey
        FROM nodes
        WHERE generation_first = ?
        ORDER BY space_id
        """,
        (generation,),
    )


def _restore_completed_state(
    connection: sqlite3.Connection,
    output_dir: Path,
) -> tuple[
    int,
    dict[str, str],
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, Any],
    list[dict[str, Any]],
    set[str],
    dict[str, Any],
]:
    completed_generation = 0
    generation_summaries: list[dict[str, Any]] = []
    while True:
        candidate = output_dir / f"G{completed_generation + 1}_generation_summary.json"
        if not candidate.exists():
            break
        with candidate.open("r", encoding="utf-8") as handle:
            summary = json.load(handle)
        if summary.get("status") != "complete":
            break
        generation_summaries.append(summary)
        completed_generation += 1

    # Runs created before parent-level checkpointing still have trustworthy
    # completed-generation summaries. Backfill those parent rows so the
    # finalized release has a complete progress ledger without replaying any
    # chemistry. A partially populated ledger is not silently repaired.
    for generation in range(1, completed_generation + 1):
        parent_count = int(
            connection.execute(
                "SELECT COUNT(*) FROM nodes WHERE generation_first = ?",
                (generation - 1,),
            ).fetchone()[0]
        )
        checkpoint_count = int(
            connection.execute(
                """
                SELECT COUNT(*)
                FROM generation_parent_progress
                WHERE generation = ? AND status = 'complete'
                """,
                (generation,),
            ).fetchone()[0]
        )
        if checkpoint_count == parent_count:
            continue
        if checkpoint_count:
            raise RuntimeError(
                "Completed generation has an incomplete parent checkpoint "
                f"ledger: G{generation}, checkpoints={checkpoint_count}, "
                f"parents={parent_count}"
            )
        # Materialize the read cursor before opening the write transaction.
        # Streaming a SELECT cursor through executemany can lock a large
        # legacy WAL database against its own checkpoint-ledger migration.
        parent_rows = list(
            connection.execute(
                """
                SELECT space_id
                FROM nodes
                WHERE generation_first = ?
                ORDER BY space_id
                """,
                (generation - 1,),
            )
        )
        connection.executemany(
            """
            INSERT INTO generation_parent_progress
            VALUES (?, ?, ?, 'complete')
            """,
            (
                (generation, source_space_id, ordinal)
                for ordinal, (source_space_id,) in enumerate(parent_rows, start=1)
            ),
        )
        connection.commit()

    # Preserve the first unfinished generation. Every parent is committed as
    # an atomic checkpoint, so an interrupted run can continue without
    # discarding hours of valid enumeration.
    first_incomplete_generation = completed_generation + 1
    connection.execute(
        "DELETE FROM derivation_events WHERE generation > ?",
        (first_incomplete_generation,),
    )
    connection.execute(
        "DELETE FROM application_audit WHERE generation > ?",
        (first_incomplete_generation,),
    )
    connection.execute(
        "DELETE FROM rejection_events WHERE generation > ?",
        (first_incomplete_generation,),
    )
    connection.execute(
        "DELETE FROM generation_parent_progress WHERE generation > ?",
        (first_incomplete_generation,),
    )
    connection.execute(
        "DELETE FROM nodes WHERE generation_first > ?",
        (first_incomplete_generation,),
    )
    connection.commit()

    full_to_space: dict[str, str] = {}
    connectivity_to_g0: dict[str, set[str]] = {}
    for space_id, generation, full_key, connectivity_key, source_ids in connection.execute(
        """
        SELECT space_id, generation_first, full_inchikey, connectivity_key,
               source_node_ids
        FROM nodes
        """
    ):
        full_to_space[full_key] = space_id
        if generation == 0:
            connectivity_to_g0.setdefault(connectivity_key, set()).update(
                token for token in source_ids.split(";") if token
            )

    parent_lookup: dict[str, set[str]] = {}
    for target_space_id, source_full_key in connection.execute(
        """
        SELECT e.target_space_id, n.full_inchikey
        FROM derivation_events e
        JOIN nodes n ON e.source_space_id = n.space_id
        """
    ):
        parent_lookup.setdefault(target_space_id, set()).add(source_full_key)

    previously_activated_rules = {
        row[0]
        for row in connection.execute(
            """
            SELECT DISTINCT grammar_rule_id
            FROM application_audit
            WHERE generation <= ?
            """,
            (completed_generation,),
        )
    }
    g0_count = connection.execute(
        "SELECT COUNT(*) FROM nodes WHERE generation_first = 0"
    ).fetchone()[0]
    g0_summary = {
        "input_records": int(g0_count),
        "valid_records": int(g0_count),
        "invalid_records": 0,
        "unique_full_stereo_structures": int(g0_count),
        "restored_from_checkpoint": True,
    }
    progress_row = connection.execute(
        """
        SELECT parent_ordinal, source_space_id
        FROM generation_parent_progress
        WHERE generation = ? AND status = 'complete'
        ORDER BY parent_ordinal DESC
        LIMIT 1
        """,
        (first_incomplete_generation,),
    ).fetchone()
    progress_source = "generation_parent_progress"
    if progress_row is None:
        legacy_source_row = connection.execute(
            """
            SELECT MAX(source_space_id)
            FROM application_audit
            WHERE generation = ?
            """,
            (first_incomplete_generation,),
        ).fetchone()
        legacy_source_id = (
            legacy_source_row[0] if legacy_source_row is not None else None
        )
        progress_row = (
            (0, legacy_source_id) if legacy_source_id is not None else None
        )
        progress_source = "legacy_application_audit_lower_bound"
    partial_generation_row_count = sum(
        int(
            connection.execute(
                f"SELECT COUNT(*) FROM {table} WHERE {column} = ?",
                (first_incomplete_generation,),
            ).fetchone()[0]
        )
        for table, column in (
            ("nodes", "generation_first"),
            ("derivation_events", "generation"),
            ("application_audit", "generation"),
            ("rejection_events", "generation"),
            ("generation_parent_progress", "generation"),
        )
    )
    partial_state = {
        "generation": first_incomplete_generation,
        "has_partial_rows": partial_generation_row_count > 0,
        "resume_after_parent_ordinal": (
            int(progress_row[0]) if progress_row is not None else 0
        ),
        "resume_after_source_space_id": (
            str(progress_row[1]) if progress_row is not None else ""
        ),
        "progress_source": progress_source,
        "partial_generation_row_count": partial_generation_row_count,
    }
    return (
        completed_generation,
        full_to_space,
        connectivity_to_g0,
        parent_lookup,
        g0_summary,
        generation_summaries,
        previously_activated_rules,
        partial_state,
    )


def _restore_partial_generation_statistics(
    connection: sqlite3.Connection,
    generation: int,
) -> tuple[Counter[str], Counter[str], set[str], int]:
    counters: Counter[str] = Counter()
    application_row = connection.execute(
        """
        SELECT COUNT(*), COALESCE(SUM(raw_product_tuple_count), 0)
        FROM application_audit
        WHERE generation = ?
        """,
        (generation,),
    ).fetchone()
    counters["matched_parent_rule_pairs"] = int(application_row[0])
    counters["raw_product_tuples"] = int(application_row[1])
    counters["unique_new_nodes"] = int(
        connection.execute(
            "SELECT COUNT(*) FROM nodes WHERE generation_first = ?",
            (generation,),
        ).fetchone()[0]
    )
    event_row = connection.execute(
        """
        SELECT COUNT(*),
               COALESCE(SUM(known_g0_full_match), 0),
               COALESCE(SUM(
                 CASE
                   WHEN known_g0_full_match = 0
                    AND known_g0_connectivity_match = 1
                   THEN 1 ELSE 0
                 END
               ), 0),
               COALESCE(SUM(immediate_reverse_cycle), 0)
        FROM derivation_events
        WHERE generation = ?
        """,
        (generation,),
    ).fetchone()
    counters["accepted_derivation_events"] = int(event_row[0])
    counters["known_G0_full_recovery_events"] = int(event_row[1])
    counters["known_G0_connectivity_only_recovery_events"] = int(event_row[2])
    counters["immediate_reverse_cycle_events"] = int(event_row[3])

    rejection_counts: Counter[str] = Counter()
    for (payload,) in connection.execute(
        """
        SELECT rejection_reason_counts
        FROM application_audit
        WHERE generation = ?
        """,
        (generation,),
    ):
        try:
            rejection_counts.update(
                {
                    str(key): int(value)
                    for key, value in json.loads(payload or "{}").items()
                }
            )
        except (TypeError, ValueError, json.JSONDecodeError):
            rejection_counts["unparseable_rejection_audit"] += 1
    activated_rules = {
        str(row[0])
        for row in connection.execute(
            """
            SELECT DISTINCT grammar_rule_id
            FROM application_audit
            WHERE generation = ?
            """,
            (generation,),
        )
    }
    maximum_suffix = 0
    for (space_id,) in connection.execute(
        "SELECT space_id FROM nodes WHERE generation_first = ?",
        (generation,),
    ):
        try:
            maximum_suffix = max(
                maximum_suffix, int(str(space_id).rsplit("_", 1)[1])
            )
        except (IndexError, ValueError):
            continue
    return counters, rejection_counts, activated_rules, maximum_suffix


def generate_chemical_space(
    grammar_path: Path,
    nodes_path: Path,
    output_dir: Path,
    *,
    max_generation: int = 3,
    max_products_per_parent_rule: int = 256,
    min_source_atom_retention: float = 0.65,
    max_abs_formal_charge: int = 2,
    resume: bool = False,
) -> dict[str, Path]:
    start = time.time()
    Chem, *_rest = require_rdkit()
    from rdkit import DataStructs
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.error")

    output_dir = ensure_dir(output_dir)
    database_path = output_dir / "taxane_reaction_grammar_space.sqlite"
    resuming_existing = bool(resume and database_path.exists())
    if database_path.exists() and not resume:
        database_path.unlink()
    connection = _initialize_database(database_path)
    grammar = read_table(grammar_path)
    rules, compile_failures = _compile_rules(grammar)
    if resuming_existing:
        (
            completed_generation,
            full_to_space,
            connectivity_to_g0,
            parent_lookup,
            g0_summary,
            generation_summaries,
            previously_activated_rules,
            partial_state,
        ) = _restore_completed_state(connection, output_dir)
    else:
        full_to_space, connectivity_to_g0, g0_summary = _load_g0(
            connection, nodes_path
        )
        completed_generation = 0
        parent_lookup = {}
        generation_summaries = []
        previously_activated_rules = set()
        partial_state = {
            "generation": 1,
            "has_partial_rows": False,
            "resume_after_parent_ordinal": 0,
            "resume_after_source_space_id": "",
            "progress_source": "fresh_run",
            "partial_generation_row_count": 0,
        }

    for generation in range(completed_generation + 1, max_generation + 1):
        generation_start = time.time()
        is_partial_resume = bool(
            resuming_existing
            and generation == partial_state["generation"]
            and partial_state["has_partial_rows"]
        )
        if is_partial_resume:
            (
                counters,
                rejection_counts,
                activated_rules_this_generation,
                new_node_counter,
            ) = _restore_partial_generation_statistics(
                connection, generation
            )
        else:
            counters = Counter()
            rejection_counts = Counter()
            activated_rules_this_generation = set()
            new_node_counter = 0
        parents = list(_generation_parents(connection, generation - 1))
        if not parents:
            generation_summaries.append(
                {
                    "generation": generation,
                    "frontier_parent_count": 0,
                    "status": "no_frontier_parents",
                }
            )
            break

        resume_after_parent_ordinal = 0
        if is_partial_resume:
            resume_after_parent_ordinal = int(
                partial_state["resume_after_parent_ordinal"]
            )
            if (
                resume_after_parent_ordinal <= 0
                and partial_state["resume_after_source_space_id"]
            ):
                parent_ordinal_by_id = {
                    row[0]: index
                    for index, row in enumerate(parents, start=1)
                }
                resume_after_parent_ordinal = parent_ordinal_by_id.get(
                    partial_state["resume_after_source_space_id"], 0
                )
        parents_to_process = parents[resume_after_parent_ordinal:]
        for parent_index, (
            source_space_id,
            source_smiles,
            source_key,
        ) in enumerate(
            parents_to_process,
            start=resume_after_parent_ordinal + 1,
        ):
            source = Chem.MolFromSmiles(source_smiles)
            if source is None:
                counters["invalid_parent_smiles"] += 1
                connection.execute(
                    """
                    INSERT OR REPLACE INTO generation_parent_progress
                    VALUES (?, ?, ?, 'complete')
                    """,
                    (generation, source_space_id, parent_index),
                )
                if parent_index % 100 == 0:
                    connection.commit()
                continue
            source_fp_pattern = _bit_vector_as_int(
                Chem.PatternFingerprint(source, fpSize=2048)
            )
            source_fp = _fingerprint(source)
            source_formula = _molecule_identity(source)["formula"]
            source_parent_keys = parent_lookup.get(source_space_id, set())

            for rule in rules:
                if rule.query_fp & ~source_fp_pattern:
                    continue
                if not source.HasSubstructMatch(rule.query):
                    continue
                activated_rules_this_generation.add(rule.grammar_rule_id)
                counters["matched_parent_rule_pairs"] += 1
                site_matches = source.GetSubstructMatches(
                    rule.query,
                    uniquify=True,
                    maxMatches=max_products_per_parent_rule + 1,
                )
                try:
                    outcomes = rule.reaction.RunReactants(
                        (source,),
                        maxProducts=max_products_per_parent_rule + 1,
                    )
                except Exception as exc:
                    reason = f"reaction_execution_failed:{type(exc).__name__}"
                    rejection_counts[reason] += 1
                    connection.execute(
                        "INSERT INTO application_audit VALUES (?, ?, ?, ?, 0, 0, 0, 1, ?, ?)",
                        (
                            generation,
                            source_space_id,
                            rule.grammar_rule_id,
                            len(site_matches),
                            "reaction_execution_failed",
                            json.dumps({reason: 1}, sort_keys=True),
                        ),
                    )
                    continue

                counters["raw_product_tuples"] += len(outcomes)
                if len(outcomes) > max_products_per_parent_rule:
                    reason = "excessive_enumeration"
                    rejection_counts[reason] += len(outcomes)
                    connection.execute(
                        "INSERT INTO application_audit VALUES (?, ?, ?, ?, ?, 0, 0, ?, ?, ?)",
                        (
                            generation,
                            source_space_id,
                            rule.grammar_rule_id,
                            len(site_matches),
                            len(outcomes),
                            len(outcomes),
                            reason,
                            json.dumps({reason: len(outcomes)}, sort_keys=True),
                        ),
                    )
                    continue

                application_rejections: Counter[str] = Counter()
                accepted_target_ids: set[str] = set()
                accepted_count = 0
                for raw_index, outcome in enumerate(outcomes, start=1):
                    product_smiles_for_audit = ""
                    reason = ""
                    if len(outcome) != 1:
                        reason = "multi_product_tuple"
                    else:
                        product = outcome[0]
                        origin_retention = _origin_retention(product, source)
                        observed_changed = _observed_changed_source_atoms(
                            source, product
                        )
                        try:
                            Chem.SanitizeMol(product)
                            product = Chem.RemoveHs(product)
                            if len(Chem.GetMolFrags(product)) != 1:
                                raise ValueError("disconnected_product")
                            if any(
                                atom.GetAtomicNum() not in ALLOWED_ATOMIC_NUMBERS
                                for atom in product.GetAtoms()
                            ):
                                raise ValueError("unsupported_element")
                            if abs(Chem.GetFormalCharge(product)) > max_abs_formal_charge:
                                raise ValueError("formal_charge_out_of_range")
                            if origin_retention < min_source_atom_retention:
                                raise ValueError("insufficient_source_atom_retention")
                            _clear_atom_maps(product)
                            identity = _molecule_identity(product)
                            product_smiles_for_audit = identity["smiles"]
                            if identity["full_inchikey"] == source_key:
                                raise ValueError("identity_no_change")
                            observed_delta = formula_delta(
                                source_formula, identity["formula"]
                            )
                            if rule.expected_delta and not atom_delta_matches(
                                observed_delta, rule.expected_delta
                            ):
                                raise ValueError("element_delta_mismatch")
                        except Exception as exc:
                            reason = str(exc) or type(exc).__name__

                    if reason:
                        application_rejections[reason] += 1
                        rejection_counts[reason] += 1
                        connection.execute(
                            "INSERT INTO rejection_events VALUES (?, ?, ?, ?, ?, ?)",
                            (
                                generation,
                                source_space_id,
                                rule.grammar_rule_id,
                                raw_index,
                                reason,
                                product_smiles_for_audit,
                            ),
                        )
                        continue

                    target_key = identity["full_inchikey"]
                    target_space_id = full_to_space.get(target_key)
                    target_is_new = target_space_id is None
                    known_full = target_key in full_to_space and target_space_id.startswith("G0_")
                    connectivity_ids = connectivity_to_g0.get(
                        identity["connectivity_key"], set()
                    )
                    known_connectivity = bool(connectivity_ids)
                    if target_is_new:
                        new_node_counter += 1
                        target_space_id = f"G{generation}_{new_node_counter:08d}"
                        connection.execute(
                            """
                            INSERT INTO nodes VALUES (
                              ?, ?, '', '', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                              0, ?, ?
                            )
                            """,
                            (
                                target_space_id,
                                generation,
                                identity["smiles"],
                                target_key,
                                identity["connectivity_key"],
                                identity["formula"],
                                identity["exact_mass"],
                                identity["heavy_atom_count"],
                                identity["clogp"],
                                identity["tpsa"],
                                identity["hbd"],
                                identity["hba"],
                                identity["rotatable_bonds"],
                                identity["ring_count"],
                                identity["fraction_csp3"],
                                identity["formal_charge"],
                                int(known_connectivity),
                                ";".join(sorted(connectivity_ids)),
                            ),
                        )
                        full_to_space[target_key] = target_space_id
                        counters["unique_new_nodes"] += 1

                    target_generation = connection.execute(
                        "SELECT generation_first FROM nodes WHERE space_id = ?",
                        (target_space_id,),
                    ).fetchone()[0]
                    target_mol = Chem.MolFromSmiles(identity["smiles"])
                    tanimoto = DataStructs.TanimotoSimilarity(
                        source_fp, _fingerprint(target_mol)
                    )
                    immediate_reverse = (
                        target_key in source_parent_keys
                        or target_space_id == source_space_id
                    )
                    connection.execute(
                        """
                        INSERT INTO derivation_events VALUES (
                          NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                        )
                        """,
                        (
                            generation,
                            source_space_id,
                            target_space_id,
                            int(target_is_new),
                            target_generation,
                            rule.grammar_rule_id,
                            rule.smarts_rule_id,
                            rule.semantic_group_id,
                            rule.reaction_type,
                            rule.evidence_layer,
                            rule.final_rule_confidence,
                            rule.expected_delta_text,
                            json.dumps(observed_delta, sort_keys=True),
                            round(origin_retention, 6),
                            observed_changed,
                            round(float(tanimoto), 6),
                            raw_index,
                            int(known_full),
                            int(known_connectivity),
                            int(immediate_reverse),
                        ),
                    )
                    accepted_count += 1
                    accepted_target_ids.add(target_space_id)
                    counters["accepted_derivation_events"] += 1
                    if known_full:
                        counters["known_G0_full_recovery_events"] += 1
                    elif known_connectivity:
                        counters["known_G0_connectivity_only_recovery_events"] += 1
                    if immediate_reverse:
                        counters["immediate_reverse_cycle_events"] += 1
                    parent_lookup.setdefault(target_space_id, set()).add(source_key)

                connection.execute(
                    "INSERT INTO application_audit VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        generation,
                        source_space_id,
                        rule.grammar_rule_id,
                        len(site_matches),
                        len(outcomes),
                        accepted_count,
                        len(accepted_target_ids),
                        sum(application_rejections.values()),
                        "processed",
                        json.dumps(application_rejections, sort_keys=True),
                    ),
                )

            connection.execute(
                """
                INSERT OR REPLACE INTO generation_parent_progress
                VALUES (?, ?, ?, 'complete')
                """,
                (generation, source_space_id, parent_index),
            )
            if parent_index % 100 == 0:
                connection.commit()
                print(
                    f"[generate-space] G{generation} parents={parent_index}/{len(parents)} "
                    f"new={counters['unique_new_nodes']} "
                    f"events={counters['accepted_derivation_events']}",
                    flush=True,
                )
        connection.commit()

        generation_summary = {
            "generation": generation,
            "frontier_parent_count": len(parents),
            "compiled_grammar_rules": len(rules),
            "activated_rules": len(activated_rules_this_generation),
            "newly_activated_rules": len(
                activated_rules_this_generation - previously_activated_rules
            ),
            "cumulative_activated_rules": len(
                activated_rules_this_generation | previously_activated_rules
            ),
            **{key: int(value) for key, value in sorted(counters.items())},
            "rejection_reason_counts": {
                key: int(value) for key, value in sorted(rejection_counts.items())
            },
            "elapsed_seconds": round(time.time() - generation_start, 3),
            "partial_generation_resumed": is_partial_resume,
            "parents_already_complete_at_start": resume_after_parent_ordinal,
            "checkpoint_progress_source": (
                partial_state["progress_source"]
                if is_partial_resume
                else "fresh_generation"
            ),
            "status": "complete",
        }
        generation_summaries.append(generation_summary)
        previously_activated_rules.update(activated_rules_this_generation)
        write_json(
            generation_summary,
            output_dir / f"G{generation}_generation_summary.json",
        )

    outputs = {
        "database": database_path,
        "nodes": output_dir / "chemical_space_nodes.tsv",
        "events": output_dir / "derivation_events.tsv",
        "application_audit": output_dir / "rule_application_audit.tsv",
        "rejections": output_dir / "rejection_events.tsv",
        "parent_progress": output_dir / "generation_parent_progress.tsv",
        "compile_failures": output_dir / "grammar_compile_failures.tsv",
        "summary": output_dir / "chemical_space_build_summary.json",
    }
    _export_query(connection, "SELECT * FROM nodes ORDER BY generation_first, space_id", outputs["nodes"])
    _export_query(connection, "SELECT * FROM derivation_events ORDER BY event_id", outputs["events"])
    _export_query(
        connection,
        "SELECT * FROM application_audit ORDER BY generation, source_space_id, grammar_rule_id",
        outputs["application_audit"],
    )
    _export_query(
        connection,
        "SELECT * FROM rejection_events ORDER BY generation, source_space_id, grammar_rule_id, raw_product_index",
        outputs["rejections"],
    )
    _export_query(
        connection,
        """
        SELECT * FROM generation_parent_progress
        ORDER BY generation, parent_ordinal
        """,
        outputs["parent_progress"],
    )
    pd.DataFrame(compile_failures).to_csv(
        outputs["compile_failures"], sep="\t", index=False
    )
    node_count = connection.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    event_count = connection.execute(
        "SELECT COUNT(*) FROM derivation_events"
    ).fetchone()[0]
    application_count = connection.execute(
        "SELECT COUNT(*) FROM application_audit"
    ).fetchone()[0]
    rejection_count = connection.execute(
        "SELECT COUNT(*) FROM rejection_events"
    ).fetchone()[0]
    parent_progress_count = connection.execute(
        "SELECT COUNT(*) FROM generation_parent_progress"
    ).fetchone()[0]
    summary = {
        "mode": "strict_iterative_reaction_grammar_derivation",
        "grammar_input": str(grammar_path),
        "nodes_input": str(nodes_path),
        "identity_policy": {
            "primary": "full_stereochemistry_aware_InChIKey",
            "auxiliary": "first_14_characters_connectivity_InChIKey",
            "connectivity_key_is_never_used_to_deduplicate_stereoisomers": True,
        },
        "generation_policy": (
            "G1 derives from unique G0 structures; each later generation derives "
            "only from the immediately preceding novel frontier."
        ),
        "compiled_grammar_rules": len(rules),
        "grammar_compile_failures": len(compile_failures),
        "G0": g0_summary,
        "generations": generation_summaries,
        "final_node_count": int(node_count),
        "final_derivation_event_count": int(event_count),
        "rule_application_audit_rows": int(application_count),
        "rejection_event_count": int(rejection_count),
        "completed_parent_checkpoints": int(parent_progress_count),
        "max_generation": max_generation,
        "resumed_from_existing_checkpoint": resuming_existing,
        "completed_generation_at_start": completed_generation,
        "partial_generation_state_at_start": partial_state,
        "max_products_per_parent_rule": max_products_per_parent_rule,
        "min_source_atom_retention": min_source_atom_retention,
        "rdkit_console_errors_suppressed": True,
        "elapsed_seconds": round(time.time() - start, 3),
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    write_json(summary, outputs["summary"])
    connection.close()
    return outputs

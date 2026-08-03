from __future__ import annotations

import json
import re
import sqlite3
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd

from .schema import FAMILY_EVIDENCE_COLUMNS, NORMALIZED_COLUMNS
from .ec import normalize_ec_numbers as normalize_ec_field, normalize_ec_number, normalize_ec_prefix
from .source_layers import infer_evidence_layer_from_record, normalize_evidence_layer
from .utils import clean_text, join_values, rel_or_abs, read_yaml, sha256_text


def _empty_norm() -> pd.DataFrame:
    return pd.DataFrame(columns=NORMALIZED_COLUMNS)


EC_RE = re.compile(r"(?<![\w.])(?:\d{1,2}|-)(?:\.(?:\d{1,3}|-|n)){3}(?![\w.])", re.IGNORECASE)


def normalize_ec(x: Any) -> str:
    vals = EC_RE.findall(clean_text(x))
    return join_values([v.replace("n", "-").replace("N", "-") for v in vals])


def normalize_rhea_ids(*texts: Any, strict: bool = True, allow_bare: bool = False) -> str:
    """Normalize Rhea identifiers without stealing numeric IDs from other databases.

    By default only IDs explicitly prefixed with RHEA are accepted. Bare 4--6 digit
    numbers are accepted only when the caller is reading a field whose schema is
    known to contain Rhea IDs, such as a Rhea database_reaction_id column or a
    RetroRules `rhea` column. This prevents values such as KEGG:R00623 from being
    converted incorrectly to RHEA:00623.
    """
    vals = []
    for text in texts:
        s = clean_text(text)
        if not s:
            continue
        vals.extend([f"RHEA:{m}" for m in re.findall(r"(?i)\bRHEA[:_\-]?(\d{4,6})\b", s)])
        if allow_bare and not strict:
            # Accept bare IDs only as standalone semicolon/comma/pipe/space-separated
            # tokens, not as the numeric suffix of KEGG:R00623, MNXR10001, etc.
            for token in re.split(r"[;,|\s]+", s):
                token = clean_text(token)
                if re.fullmatch(r"\d{4,6}", token):
                    vals.append(f"RHEA:{token}")
    return join_values(vals)


def _row_hash(row: dict[str, Any]) -> str:
    return sha256_text(json.dumps(row, ensure_ascii=False, sort_keys=True, default=str))[:20]


def _base_out(dataset: dict, source_file: str, parser_name: str, row_i: int | str = "") -> dict[str, str]:
    source_database = clean_text(dataset.get("source_database", "Unknown")) or "Unknown"
    explicit_layer = clean_text(dataset.get("evidence_layer", ""))
    out = {c: "" for c in NORMALIZED_COLUMNS}
    out["source_database"] = source_database
    out["evidence_layer"] = explicit_layer or "Unknown"
    out["source_file"] = source_file
    out["parser_name"] = parser_name
    out["raw_row_index"] = clean_text(row_i)
    return out


def _finalize_row(out: dict[str, Any], raw: dict[str, Any] | None = None) -> dict[str, str]:
    if not clean_text(out.get("reaction_smiles", "")) and clean_text(out.get("substrate_smiles", "")) and clean_text(out.get("product_smiles", "")):
        out["reaction_smiles"] = f"{out['substrate_smiles']}>>{out['product_smiles']}"
    out["ec_numbers"] = normalize_ec(out.get("ec_numbers", ""))
    out["template_ec_candidates"] = normalize_ec(out.get("template_ec_candidates", ""))
    out["database_ec_candidates"] = normalize_ec(out.get("database_ec_candidates", ""))
    out["ec_prior_candidates"] = normalize_ec(out.get("ec_prior_candidates", ""))
    out["rhea_ids"] = normalize_rhea_ids(out.get("rhea_ids", ""), strict=True)
    out["evidence_layer"] = normalize_evidence_layer(infer_evidence_layer_from_record(out, out.get("evidence_layer", "Unknown")))
    if not clean_text(out.get("record_id", "")):
        rid = clean_text(out.get("source_reaction_id", "")) or clean_text(out.get("raw_row_index", ""))
        out["record_id"] = f"{out.get('source_database', 'Unknown')}:{rid}"
    out["source_row_hash"] = _row_hash(raw or out)
    return {c: clean_text(out.get(c, "")) for c in NORMALIZED_COLUMNS}


def _first_present(row: pd.Series, names: list[str]) -> str:
    for n in names:
        if n in row.index and clean_text(row.get(n, "")):
            return clean_text(row.get(n, ""))
    return ""


def _apply_column_map(df: pd.DataFrame, dataset: dict, source_file: str, parser_name: str) -> pd.DataFrame:
    cmap = dataset.get("column_map", {}) or {}
    rows = []
    for i, row in df.fillna("").iterrows():
        out = _base_out(dataset, source_file, parser_name, i)
        raw = row.to_dict()
        for std_col, src_col in cmap.items():
            if src_col in df.columns and std_col in out:
                out[std_col] = clean_text(row.get(src_col, ""))
        passthrough = [
            "source_reaction_id", "reaction_smiles", "reaction_smarts", "substrate_smiles", "product_smiles",
            "reaction_equation", "ec_numbers", "template_ec_candidates", "database_ec_candidates", "ec_prior_candidates",
            "direction", "is_reversible", "reaction_type_source", "reaction_subtype_source", "cofactor_or_donor_class",
            "source_evidence_text", "protein_ids", "enzyme_name", "rhea_ids", "kegg_ids", "metanetx_ids",
        ]
        for std_col in passthrough:
            if not out.get(std_col) and std_col in df.columns:
                out[std_col] = clean_text(row.get(std_col, ""))
        rows.append(_finalize_row(out, raw))
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS)


def load_taxol_pathway_csv(path: str | Path, dataset: dict) -> pd.DataFrame:
    """Load curated known Taxol pathway reactions as source records.

    They are exact source/anchor records. They become generalized predictive rules only if the input row also supplies
    a valid reaction SMARTS or an upstream tool later derives a validated template.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    rows = []
    for i, row in df.fillna("").iterrows():
        out = _base_out({**dataset, "source_database": dataset.get("source_database", "TaxolKnownPathway_Curated"), "evidence_layer": dataset.get("evidence_layer", "T1_Bio_Core")}, str(path), "taxol_pathway_csv", i)
        enzyme = clean_text(row.get("Enzyme", row.get("enzyme", "")))
        out["source_reaction_id"] = enzyme or f"taxol_pathway_{i+1}"
        out["curated_taxol_anchor"] = "true"
        out["curated_pathway_name"] = clean_text(dataset.get("curated_pathway_name", "")) or "TaxolKnownPathway"
        out["curated_pathway_step_id"] = clean_text(row.get("step_id", row.get("Step", row.get("step", "")))) or f"taxol_pathway_{i+1}"
        out["enzyme_name"] = enzyme
        out["substrate_smiles"] = clean_text(row.get("Substrate", row.get("substrate_smiles", "")))
        out["product_smiles"] = clean_text(row.get("Product", row.get("product_smiles", "")))
        out["reaction_smarts"] = clean_text(row.get("reaction_smarts", row.get("ReactionSMARTS", "")))
        out["protein_sequence"] = clean_text(row.get("Protein sequence", row.get("protein_sequence", "")))
        out["ec_numbers"] = normalize_ec(row.get("EC", row.get("ec_numbers", "")))
        out["direction"] = clean_text(row.get("direction", "forward")) or "forward"
        out["is_reversible"] = clean_text(row.get("is_reversible", "false")) or "false"
        out["source_evidence_text"] = enzyme
        rows.append(_finalize_row(out, row.to_dict()))
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS)


def load_generic_reaction_table(path: str | Path, dataset: dict) -> pd.DataFrame:
    delimiter = dataset.get("delimiter") or ("\t" if str(path).lower().endswith((".tsv", ".txt")) else ",")
    max_records = dataset.get("max_records")
    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False, low_memory=False, nrows=int(max_records) if max_records else None)
    return _apply_column_map(df, dataset, str(path), dataset.get("parser", "generic_reaction_table"))


def load_reaction_smiles_lines(path: str | Path, dataset: dict) -> pd.DataFrame:
    rows = []
    max_records = int(dataset.get("max_records") or 0)
    source_dataset = clean_text(dataset.get("source_dataset", "")) or Path(path).stem
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if max_records and len(rows) >= max_records:
                break
            rxn = clean_text(line)
            if not rxn or ">>" not in rxn:
                continue
            out = _base_out(dataset, str(path), "reaction_smiles_lines", line_no)
            out["source_reaction_id"] = f"{source_dataset}:line_{line_no}"
            out["reaction_smiles"] = rxn
            out["source_evidence_text"] = join_values([
                clean_text(dataset.get("source_evidence_text", "")),
                "raw_reaction_smiles_line",
                source_dataset,
            ])
            rows.append(_finalize_row(out, {"reaction_smiles": rxn, "line_no": line_no, "source_dataset": source_dataset}))
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS) if rows else _empty_norm()


def load_reaction_smiles_zip(path: str | Path, dataset: dict) -> pd.DataFrame:
    rows = []
    max_records = int(dataset.get("max_records") or 0)
    members = dataset.get("members") or dataset.get("zip_members") or []
    if isinstance(members, str):
        members = [m.strip() for m in re.split(r"[;,]", members) if m.strip()]
    with zipfile.ZipFile(path) as zf:
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if members:
            names = [n for n in names if n in set(members)]
        for name in names:
            source_dataset = clean_text(dataset.get("source_dataset", "")) or Path(name).stem
            with zf.open(name) as handle:
                for line_no, raw in enumerate(handle, start=1):
                    if max_records and len(rows) >= max_records:
                        break
                    rxn = clean_text(raw.decode("utf-8", errors="replace"))
                    if not rxn or ">>" not in rxn:
                        continue
                    out = _base_out(dataset, f"{path}:{name}", "reaction_smiles_zip", line_no)
                    out["source_reaction_id"] = f"{source_dataset}:{name}:line_{line_no}"
                    out["reaction_smiles"] = rxn
                    out["source_evidence_text"] = join_values([
                        clean_text(dataset.get("source_evidence_text", "")),
                        "raw_reaction_smiles_zip",
                        name,
                    ])
                    rows.append(_finalize_row(out, {"reaction_smiles": rxn, "line_no": line_no, "zip_member": name}))
                if max_records and len(rows) >= max_records:
                    break
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS) if rows else _empty_norm()


def load_reaction_evidence_csv(path: str | Path, dataset: dict) -> pd.DataFrame:
    max_records = dataset.get("max_records")
    df = pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False, nrows=int(max_records) if max_records else None)
    rows = []
    for i, row in df.fillna("").iterrows():
        row_db = clean_text(row.get("source_database", "")) or clean_text(dataset.get("source_database", "")) or "Unknown"
        out = _base_out({**dataset, "source_database": row_db}, str(path), "reaction_evidence_csv", i)
        out["source_reaction_id"] = clean_text(row.get("source_reaction_id", "")) or clean_text(row.get("database_reaction_id", "")) or f"{row_db}:{i+1}"
        out["reaction_smiles"] = clean_text(row.get("reaction_smiles", ""))
        out["reaction_equation"] = clean_text(row.get("reaction_equation", ""))
        out["ec_numbers"] = clean_text(row.get("ec_numbers", ""))
        out["kegg_ids"] = clean_text(row.get("kegg_reaction_ids", ""))
        if row_db.lower() == "kegg" and not out["kegg_ids"]:
            out["kegg_ids"] = clean_text(row.get("database_reaction_id", ""))
        out["metanetx_ids"] = clean_text(row.get("mnx_reaction_ids", "")) or (clean_text(row.get("database_reaction_id", "")) if row_db.lower() in {"metanetx", "mnx"} else "")
        # Do not parse bare numbers from cross-references: KEGG:R00623 must not become RHEA:00623.
        prefixed_rhea = normalize_rhea_ids(row.get("source_reaction_id", ""), row.get("cross_references", ""), strict=True)
        bare_rhea = normalize_rhea_ids(row.get("database_reaction_id", ""), strict=False, allow_bare=(row_db.lower() == "rhea"))
        out["rhea_ids"] = join_values([prefixed_rhea, bare_rhea])
        out["source_evidence_text"] = join_values([
            clean_text(row.get("cross_references", "")),
            f"is_balanced={clean_text(row.get('is_balanced', ''))}",
            f"is_transport={clean_text(row.get('is_transport', ''))}",
            clean_text(row.get("source_file", "")),
        ])
        rows.append(_finalize_row(out, row.to_dict()))
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS) if rows else _empty_norm()


def load_mollink_stratified_table(path: str | Path, dataset: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep=dataset.get("delimiter", "\t"), dtype=str, keep_default_na=False, low_memory=False)
    rows = []
    for i, row in df.fillna("").iterrows():
        src_db = _first_present(row, ["source_database", "template_sources", "source_databases"]) or dataset.get("source_database", "MolLink")
        explicit_layer = _first_present(row, ["template_source_layer", "evidence_layer_best", "evidence_layer", "template_layer"])
        out = _base_out({**dataset, "source_database": src_db, "evidence_layer": explicit_layer or dataset.get("evidence_layer", "")}, str(path), "mollink_stratified_table", i)
        out["source_reaction_id"] = _first_present(row, ["source_reaction_ids", "source_reaction_id", "template_id", "network_rule_uid", "rule_id"])
        out["reaction_smarts"] = _first_present(row, ["reaction_smarts", "forward_reaction_smarts", "template_smarts", "retro_template"])
        out["reaction_smiles"] = _first_present(row, ["reaction_smiles", "representative_reaction_smiles", "canonical_reaction_smiles"])
        out["substrate_smiles"] = _first_present(row, ["substrate_smiles", "main_substrate_smiles", "example_substrate_smiles"])
        out["product_smiles"] = _first_present(row, ["product_smiles", "main_product_smiles", "example_product_smiles"])
        # Keep EC evidence layers separate. `candidate_ec_numbers` from an upstream table is ambiguous,
        # so it is retained as lower-strength database/cross-reference evidence unless more specific fields exist.
        out["ec_numbers"] = _first_present(row, ["source_ec_candidates", "source_ec_numbers", "primary_ec", "ec_numbers"])
        out["template_ec_candidates"] = _first_present(row, ["template_ec_candidates"])
        out["database_ec_candidates"] = join_values([
            _first_present(row, ["database_ec_candidates"]),
            _first_present(row, ["candidate_ec_numbers"]),
        ])
        out["ec_prior_candidates"] = _first_present(row, ["ec_prior_candidates"])
        out["rhea_ids"] = _first_present(row, ["rhea_ids", "rhea_reaction_ids_clean"])
        out["kegg_ids"] = _first_present(row, ["kegg_reaction_ids", "kegg_reaction_ids_clean", "kegg_ids"])
        out["metanetx_ids"] = _first_present(row, ["mnx_reaction_ids", "mnx_reaction_ids_clean", "metanetx_ids"])
        out["reaction_type_source"] = _first_present(row, ["reaction_type", "reaction_class", "reaction_type_label", "reaction_type_id"])
        out["reaction_subtype_source"] = _first_present(row, ["reaction_subtype", "reaction_subclass"])
        out["source_evidence_text"] = _first_present(row, ["template_layer_reason", "source_evidence_class", "origin_kind", "reaction_class", "notes"])
        out["direction"] = _first_present(row, ["allowed_direction", "direction"])
        out["protein_ids"] = _first_present(row, ["protein_ids", "uniprot_accessions", "uniprot_ids"])
        rows.append(_finalize_row(out, row.to_dict()))
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS)


def load_retrorules_table(path: str | Path, dataset: dict) -> pd.DataFrame:
    return load_generic_reaction_table(path, dataset)


def load_retrorules_sqlite(path: str | Path, dataset: dict) -> pd.DataFrame:
    con = sqlite3.connect(path)
    try:
        try:
            max_records = int(dataset.get("max_records") or 0)
            chunk_size = int(dataset.get("sqlite_chunk_size") or 5000)
            ids = pd.read_sql_query(
                "SELECT DISTINCT smarts_id FROM rules WHERE smarts_id IS NOT NULL ORDER BY smarts_id",
                con,
            )["smarts_id"].astype(str).tolist()
            if max_records > 0:
                ids = ids[:max_records]
            rows = []
            for start in range(0, len(ids), chunk_size):
                chunk_ids = ids[start : start + chunk_size]
                if not chunk_ids:
                    continue
                placeholders = ",".join(["?"] * len(chunk_ids))
                # Aggregate inside SQLite before materializing rows in pandas.
                # The RetroRules MVC database stores repeated SMARTS rules
                # across many reaction/substrate contexts; chunking by the
                # indexed smarts_id keeps memory bounded for large databases.
                query = f"""
                WITH ec_agg AS (
                    SELECT reaction_id, group_concat(DISTINCT ec_number) AS ec_numbers
                    FROM ec_reactions
                    GROUP BY reaction_id
                ),
                base AS (
                    SELECT
                        rules.smarts_id AS smarts_id,
                        smarts.smarts_string AS reaction_smarts,
                        rules.reaction_id AS reaction_id,
                        rules.substrate_id AS substrate_id,
                        rules.diameter AS diameter,
                        rules.isStereo AS is_stereo,
                        rules.direction AS direction,
                        rules.score AS retrorules_score,
                        reactions.mnxr AS mnx_reaction_id,
                        reactions.kegg AS kegg_reaction_id,
                        reactions.metacyc AS metacyc_reaction_id,
                        reactions.rhea AS rhea_reaction_id,
                        reactions.reactome AS reactome_reaction_id,
                        reactions.brenda AS brenda_reaction_id,
                        ec_agg.ec_numbers AS ec_numbers
                    FROM rules
                    JOIN smarts ON rules.smarts_id = smarts.id
                    LEFT JOIN reactions ON rules.reaction_id = reactions.id
                    LEFT JOIN ec_agg ON rules.reaction_id = ec_agg.reaction_id
                    WHERE rules.smarts_id IN ({placeholders})
                )
                SELECT
                    smarts_id,
                    max(reaction_smarts) AS reaction_smarts,
                    group_concat(DISTINCT reaction_id) AS reaction_id,
                    group_concat(DISTINCT substrate_id) AS substrate_id,
                    group_concat(DISTINCT diameter) AS diameter,
                    group_concat(DISTINCT is_stereo) AS is_stereo,
                    group_concat(DISTINCT direction) AS direction,
                    max(retrorules_score) AS retrorules_score,
                    group_concat(DISTINCT mnx_reaction_id) AS mnx_reaction_id,
                    group_concat(DISTINCT kegg_reaction_id) AS kegg_reaction_id,
                    group_concat(DISTINCT metacyc_reaction_id) AS metacyc_reaction_id,
                    group_concat(DISTINCT rhea_reaction_id) AS rhea_reaction_id,
                    group_concat(DISTINCT reactome_reaction_id) AS reactome_reaction_id,
                    group_concat(DISTINCT brenda_reaction_id) AS brenda_reaction_id,
                    group_concat(DISTINCT ec_numbers) AS ec_numbers
                FROM base
                GROUP BY smarts_id
                ORDER BY smarts_id
                """
                df = pd.read_sql_query(query, con, params=chunk_ids).fillna("")
                for i, row in df.fillna("").iterrows():
                    smarts = clean_text(row.get("reaction_smarts", ""))
                    if not smarts or ">>" not in smarts:
                        continue
                    out = _base_out({**dataset, "source_database": dataset.get("source_database", "RetroRules")}, str(path), "retrorules_sqlite", i)
                    ecs = normalize_ec_field(row.get("ec_numbers", ""))
                    out["source_reaction_id"] = "RetroRules:" + clean_text(row.get("reaction_id", ""))
                    out["reaction_smarts"] = smarts
                    out["ec_numbers"] = ecs
                    out["template_ec_candidates"] = ecs
                    out["direction"] = clean_text(row.get("direction", ""))
                    out["rhea_ids"] = normalize_rhea_ids(row.get("rhea_reaction_id", ""), strict=False, allow_bare=True)
                    out["kegg_ids"] = clean_text(row.get("kegg_reaction_id", ""))
                    out["metanetx_ids"] = clean_text(row.get("mnx_reaction_id", ""))
                    out["source_evidence_text"] = join_values([
                        "retrorules_sqlite_join",
                        f"smarts_id={clean_text(row.get('smarts_id', ''))}",
                        f"diameter={clean_text(row.get('diameter', ''))}",
                        f"is_stereo={clean_text(row.get('is_stereo', ''))}",
                        f"score={clean_text(row.get('retrorules_score', ''))}",
                        f"metacyc={clean_text(row.get('metacyc_reaction_id', ''))}",
                        f"reactome={clean_text(row.get('reactome_reaction_id', ''))}",
                        f"brenda={clean_text(row.get('brenda_reaction_id', ''))}",
                    ])
                    hash_payload = {
                        "source_database": out.get("source_database", "RetroRules"),
                        "smarts_id": clean_text(row.get("smarts_id", "")),
                        "reaction_id": clean_text(row.get("reaction_id", "")),
                        "substrate_id": clean_text(row.get("substrate_id", "")),
                        "diameter": clean_text(row.get("diameter", "")),
                        "direction": clean_text(row.get("direction", "")),
                    }
                    rows.append(_finalize_row(out, hash_payload))
            if rows:
                return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS)
        except Exception:
            # Fall through to schema-agnostic loading for non-standard RetroRules exports.
            pass

        tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", con)["name"].tolist()
        frames = []
        for t in tables:
            try:
                df = pd.read_sql_query(f"SELECT * FROM {t}", con)
            except Exception:
                continue
            low = {c.lower(): c for c in df.columns}
            if not any(k in low for k in ["reaction_smarts", "smarts", "rule_smarts", "template"]):
                continue
            cmap = dict(dataset.get("column_map", {}) or {})
            for cand in ["reaction_smarts", "smarts", "rule_smarts", "template"]:
                if cand in low and "reaction_smarts" not in cmap:
                    cmap["reaction_smarts"] = low[cand]
                    break
            for cand in ["rule_id", "id", "reaction_id", "reaction"]:
                if cand in low and "source_reaction_id" not in cmap:
                    cmap["source_reaction_id"] = low[cand]
                    break
            for cand in ["ec", "ec_number", "ec_numbers"]:
                if cand in low and "ec_numbers" not in cmap:
                    cmap["ec_numbers"] = low[cand]
                    break
            for cand in ["diameter", "radius"]:
                if cand in low and "source_evidence_text" not in cmap:
                    cmap["source_evidence_text"] = low[cand]
                    break
            ds = dict(dataset)
            ds["column_map"] = cmap
            frames.append(_apply_column_map(df, ds, f"{path}:{t}", "retrorules_sqlite"))
        return pd.concat(frames, ignore_index=True) if frames else _empty_norm()
    finally:
        con.close()


def load_datasets(datasets_yaml: str | Path, override_taxol_path: str | Path | None = None) -> tuple[pd.DataFrame, list[str]]:
    cfg = read_yaml(datasets_yaml)
    base = Path(datasets_yaml).resolve().parent
    frames, used = [], []
    for name, ds0 in (cfg.get("datasets", {}) or {}).items():
        ds = dict(ds0 or {})
        if not ds.get("enabled", False):
            continue
        parser = ds.get("parser", "generic_reaction_table")
        path0 = override_taxol_path if parser == "taxol_pathway_csv" and override_taxol_path else ds.get("path")
        if not path0:
            continue
        path = rel_or_abs(path0, base)
        if not path.exists():
            raise FileNotFoundError(f"Dataset {name} path not found: {path}")
        used.append(str(path))
        if parser == "taxol_pathway_csv":
            frames.append(load_taxol_pathway_csv(path, ds))
        elif parser in {"generic_reaction_table", "generic_template_table"}:
            frames.append(load_generic_reaction_table(path, ds))
        elif parser in {"reaction_smiles_lines", "bionavi_reaction_lines"}:
            frames.append(load_reaction_smiles_lines(path, ds))
        elif parser in {"reaction_smiles_zip", "bionavi_reaction_zip"}:
            frames.append(load_reaction_smiles_zip(path, ds))
        elif parser in {"reaction_evidence_csv", "database_reaction_evidence_csv"}:
            frames.append(load_reaction_evidence_csv(path, ds))
        elif parser == "retrorules_table":
            frames.append(load_retrorules_table(path, ds))
        elif parser == "retrorules_sqlite":
            frames.append(load_retrorules_sqlite(path, ds))
        elif parser in {"mollink_stratified_table", "mollink_network_ready_table", "mollink_stratified", "mollink_network_ready"}:
            frames.append(load_mollink_stratified_table(path, ds))
        else:
            raise ValueError(f"Unknown parser for dataset {name}: {parser}")
    return (pd.concat(frames, ignore_index=True).fillna("") if frames else _empty_norm(), used)


def family_evidence_paths(datasets_yaml: str | Path) -> list[str]:
    cfg = read_yaml(datasets_yaml)
    fe = cfg.get("family_evidence", {}) or {}
    if not fe.get("enabled", False) or not fe.get("path"):
        return []
    return [str(rel_or_abs(fe.get("path"), Path(datasets_yaml).resolve().parent))]


def load_family_evidence(datasets_yaml: str | Path) -> pd.DataFrame:
    cfg = read_yaml(datasets_yaml)
    fe = cfg.get("family_evidence", {}) or {}
    if not fe.get("enabled", False) or not fe.get("path"):
        return pd.DataFrame(columns=FAMILY_EVIDENCE_COLUMNS)
    path = rel_or_abs(fe.get("path"), Path(datasets_yaml).resolve().parent)
    df = pd.read_csv(path, sep=fe.get("delimiter", "\t"), dtype=str, keep_default_na=False).fillna("")
    # Allow old schema: family -> primary_family.
    if "primary_family" not in df.columns and "family" in df.columns:
        df["primary_family"] = df["family"]
    if "family" not in df.columns and "primary_family" in df.columns:
        df["family"] = df["primary_family"]
    if "match_value" not in df.columns:
        if "source_reaction_id" in df.columns:
            df["match_value"] = df["source_reaction_id"]
            df["match_type"] = df.get("match_type", "source_reaction_id")
        elif "ec_number" in df.columns:
            df["match_value"] = df["ec_number"]
            df["match_type"] = df.get("match_type", "ec_prefix")
    if "evidence_id" not in df.columns:
        df["evidence_id"] = [f"FAM_EVID_{i+1:09d}" for i in range(len(df))]
    for c in FAMILY_EVIDENCE_COLUMNS:
        if c not in df.columns:
            df[c] = ""
    # Normalize EC-like columns and EC match values so family evidence joins are deterministic.
    if "ec_number" in df.columns:
        df["ec_number"] = df["ec_number"].map(normalize_ec_field)
    if "match_type" in df.columns and "match_value" in df.columns:
        mt = df["match_type"].astype(str).str.lower()
        exact_mask = mt.isin(["ec", "ec_exact"])
        prefix_mask = mt.eq("ec_prefix")
        df.loc[exact_mask, "match_value"] = df.loc[exact_mask, "match_value"].map(normalize_ec_number)
        df.loc[prefix_mask, "match_value"] = df.loc[prefix_mask, "match_value"].map(normalize_ec_prefix)
    return df[FAMILY_EVIDENCE_COLUMNS]

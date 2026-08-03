# Publication workflows

The scripts under `publication_workflows/` are presentation-only copies used
for the final article release. They read frozen source tables and do not build
reaction rules or enumerate chemical space.

After extracting the companion data record, pass its release root explicitly
to each workflow. For example:

```bash
python publication_workflows/main_figures/Figure_1/build_figure1_redesign_v6.py \
  --work /path/to/extracted/article_release \
  --output reproduced_figures/Figure_1

python publication_workflows/main_figures/Figure_2/build_figure2_v12.py \
  --work /path/to/extracted/article_release \
  --output reproduced_figures/Figure_2

python publication_workflows/main_figures/Figure_3/build_figure3_v7.py \
  --work /path/to/extracted/article_release \
  --output reproduced_figures/Figure_3

python publication_workflows/main_figures/Figure_4/build_figure4_redesigned_v8.py \
  --source-data /path/to/extracted/legacy_figure_source_data \
  --output reproduced_figures/Figure_4
```

Each workflow writes its panel-level source tables, numerical audit, and
vector/raster exports. The released source-data manifest identifies the exact
input files for every final panel.

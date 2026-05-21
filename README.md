# flexpipe-InfluenzaA

Nextstrain pipeline for genomic epidemiology of Influenza A virus (H1N1 and H3N2), segments HA and NA. This tool is derived from [flexpipe](https://github.com/InstitutoTodosPelaSaude/flexpipe) and adapted for Influenza A, similarly to how flexpipe-RSV was created for respiratory syncytial virus.

This repository contains all essential files to generate Influenza A Nextstrain builds for the HA and NA segments of both H1N1 and H3N2 subtypes. Using this pipeline, users can perform genomic epidemiology analyses, visualize phylogeographic results, and track Influenza A spread based on genomic data and associated metadata.

## Getting Started

To run this pipeline, see the instructions available in the original [flexpipe repository](https://github.com/InstitutoTodosPelaSaude/flexpipe), which covers Unix CLI navigation, installation of a Nextstrain environment with conda/mamba, and a step-by-step tutorial on generating a Nextstrain build (preparing, aligning, and visualizing genomic data).

---

## Builds

| Build | Segment | Subtype | Reference | Nextclade dataset |
|-------|---------|---------|-----------|-------------------|
| H1N1_HA | HA | H1N1 | CY121680 | `flu-h1n1-ha-CY121680` |
| H1N1_NA | NA | H1N1 | — | `flu-h1n1-na` |
| H3N2_HA | HA | H3N2 | NC_007366 | `flu-h3n2-ha-EPI1857216` |
| H3N2_NA | NA | H3N2 | — | `flu-h3n2-na` |

Each build lives in its own subdirectory with independent `config/`, `data/`, `ingest/`, `phylogenetic/`, and `scripts/` folders.

---

## Pipeline Overview

```
fetch_ncbi
    └── merge_local_sequences  (ITpS sequences + NCBI)
            └── viralqc        (Nextclade QC + clade assignment)
                    └── curate_qc  (normalisation, dedup, filters)
                            └── prepare  (subsampling)
                                    ├── coordinates  (geocoding → latlongs.tsv)
                                    ├── generate_name2hue  (colour palette)
                                    └── colours  (colour_scheme.tsv)
                                            └── [phylogenetic/Snakefile]
                                                    align → mask → tree → refine
                                                    → ancestral → translate → traits
                                                    → clades → export → auspice/results.json
```

---

## Stage 1 — Ingest

Sequences and metadata are fetched from **NCBI Virus** using Entrez, filtered by taxonomy ID and segment-specific search terms.

| Parameter | H1N1 HA | H1N1 NA | H3N2 HA | H3N2 NA |
|-----------|---------|---------|---------|---------|
| `ncbi.taxid` | 11320 | 11320 | 11320 | 11320 |
| `ncbi.genome_size` (bp) | 1752 | 1434 | 1762 | 1467 |
| `ncbi.min_length` | 80% | 80% | 80% | 80% |
| `ncbi.max_length` | 110% | 110% | 110% | 110% |
| `ncbi.min_date` | 2015-01-01 | 2015-01-01 | 2015-01-01 | 2015-01-01 |

Local ITpS sequences (in `data/new_sequences.fasta` + `data/new_metadata.tsv`) are merged with NCBI sequences at this stage.

---

## Stage 2 — QC and Curation

**Nextclade** (via ViralQC wrapper) assigns clades and quality scores to every sequence. The `curate.py` script then:

- Renames and standardises metadata fields (`strain`, `date`, `country`, `division`, `location`, `data_use`, `clade`)
- Truncates clade names to a configurable number of hierarchy levels (`clade_levels`) for display grouping
- Normalises the `host` field to canonical names (`human`, `swine`, `avian`, `ferret`, `dog`, etc.)
- Infers geographic regions from country names
- Deduplicates sequences, preferring local ITpS records

**QC filters** applied by `augur filter`:

| Parameter | Value |
|-----------|-------|
| `qc.nextclade_status` | `good`, `mediocre` |
| `qc.min_coverage` | 0.80 |
| Required columns | `strain`, `date`, `country`, `clade` |

Sequences with Nextclade status `bad` or coverage below 80% are discarded.

### Clade truncation

| Build | `clade_levels` | Example |
|-------|---------------|---------|
| H1N1 HA | 2 | `D.3.1.1` → `D.3` |
| H1N1 NA | 2 | `6B.1A.5a.2a` → `6B.1A` |
| H3N2 HA | 1 | `J.2.1.1` → `J` |
| H3N2 NA | 2 | `N2` clades → first 2 levels |

---

## Stage 3 — Subsampling

Controlled by `config/subsample.yaml`. Strategy: **focal** sequences (ITpS + all Brazil) are always kept in full; **context** sequences (outside Brazil) are subsampled by country × year × clade_truncated.

```yaml
defaults:
  min_date: 2015

samples:
  focal:
    query: "(source == 'ITpS') or (country == 'Brazil')"

  context:
    group_by: [country, year, clade_truncated]
    sequences_per_group: 5
    exclude_where:
      - "country=Brazil"
      - "clade_truncated="
      - "date="
```

Adjust `sequences_per_group` and `min_date` to control dataset size and temporal depth.

---

## Stage 4 — Coordinates and Colours

**Coordinates**: `get_coordinates.py` queries Nominatim (OpenStreetMap) to geocode `country`, `division`, and `location` fields. Results are cached in `config/cache_coordinates.tsv` to avoid redundant API calls. The output `config/latlongs.tsv` is consumed by `augur export`.

**Colours**: `generate_name2hue.py` assigns hues to each unique value in `clade_truncated`, `host`, `source`, `data_use`, and geographic columns. `colour_maker.py` produces the final `config/colour_scheme.tsv`.

Colour columns configured in `config.yaml`:

```yaml
colours:
  clade: "clade_truncated clade"
  geo:   "region country division location"
  host:  "host"
  source: "source"
  data_use: "data_use"
```

---

## Stage 5 — Phylogenetic

Run separately after ingest completes:

```bash
snakemake --snakefile phylogenetic/Snakefile --cores 8
```

Steps: `align` (MAFFT) → `mask` → `tree` (IQ-TREE 3 UFBoot) → `refine` (TreeTime) → `ancestral` → `translate` → `traits` → `clades` → `export` → `auspice/results.json`

### Key phylogenetic parameters (`config.yaml`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `parameters.model` | `MFP` | ModelFinder Plus — auto-selects best substitution model |
| `parameters.ufboot` | `1000` | Ultrafast bootstrap replicates |
| `parameters.root` | `least-squares` | Root method for time-calibrated tree |
| `parameters.coalescent` | `skyline` | Effective population size model in TreeTime |
| `parameters.date_inference` | `marginal` | Marginal date inference for ambiguous dates |
| `parameters.divergence_units` | `mutations` | Branch length units in timetree |
| `parameters.clock_filter_iqd` | `4` | IQD filter for clock outliers |
| `parameters.ancestral_inference` | `joint` | Joint ancestral reconstruction |
| `parameters.mask_5prime` | `1` | Bases masked at 5′ end |
| `parameters.mask_3prime` | `1` | Bases masked at 3′ end |
| `options.threads` | `8` | Threads for MAFFT and IQ-TREE |
| `traits.columns` | `country division location clade` | Columns for ancestral trait inference |

---

## Configuration Files

| File | Purpose |
|------|---------|
| `config/config.yaml` | All pipeline parameters (NCBI, QC, phylogenetic, colours) |
| `config/subsample.yaml` | Subsampling strategy and group sizes |
| `config/auspice_config.json` | Auspice display settings (colorings, filters, panels) |
| `config/reference.gb` | Segment reference sequence in GenBank format |
| `config/clades.tsv` | Clade definitions for `augur clades` |
| `config/keep.txt` | Strains to always include |
| `config/ignore.txt` | Strains to always exclude |
| `data/new_sequences.fasta` | Local ITpS sequences |
| `data/new_metadata.tsv` | Local ITpS metadata |

---

## Author

**Thales Bermann** — Instituto Todos pela Saúde (ITpS)
✉️ [thalesbermann@gmail.com](mailto:thalesbermann@gmail.com)

---

## License

This project is licensed under the [MIT License](LICENSE).

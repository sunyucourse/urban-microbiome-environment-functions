# Environmental correlates of urban microbiome functional potential

Processed-data reproducibility materials for **Environmental correlates of functional
potential in the global urban microbiome: cross-cohort evidence from 271 cities**.

## Scope and scientific interpretation

This repository reproduces selected **post-profiling statistical analyses** from
saved, real derived indices and environmental tables. It is **not an end-to-end
FASTQ-to-results pipeline**. Surface profiling, QC, upstream index construction,
environmental raster extraction and discovery category-PC1 construction are not
rerun. Saved category PC1 scores and factor-retention decisions are explicit inputs.

The discovery analysis uses 2,049 included surface samples from 41 cities; availability
is index-specific (not all indices were measured in all samples). Sewage city indices
represent 1,155 source samples from 230 cities. The package does not reconstruct
sewage profiles from raw reads or claim identical profiling methods across cohorts.

The central finding is selective cross-cohort consistency, not universal replication:
29/37 refers only to direction agreement for AMR depth; 185 pairs cover five
indices and 37 environmental factors. These are correlated pairs, so the reported
binomial/concordance p-values are nominal. Signed discovery latitude was historically
mapped to sewage absolute latitude; this mapping is retained to reproduce the saved
results, but it is not strict replication of the same exposure. AMR class diversity
and sewage gene richness are imperfectly matched constructs. Pathogen indices use
different species/genus definitions. LMG shares are conditional, in-sample R² allocations,
not causal effects, uniquely attributable fractions, or out-of-sample performance.

## Run

Use Python 3.11 or a compatible Python 3 environment. Exact versions tested in the
authoring environment are in `requirements.txt` and `validation/environment.json`.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead:
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python run_analysis.py --output results
```

No network downloads, API credentials, reference databases or HPC are needed.
Use a fresh output directory for each run. `python run_analysis.py --pilot --output
pilot_results` runs only the primary/sensitivity screens and factor two-stage tests.
The full run automatically compares regenerated tables against 14 historical
reference CSVs with numerical tolerances (rtol 1e-7, atol 1e-10).

## Contents

| Location | Contents |
| --- | --- |
| `data/indices/` | Saved discovery sample/city indices, including QC and source fields |
| `data/master_env_city.csv` | Discovery city environmental predictors |
| `data/city_coords_metasub.csv` | Coordinates used for spatial eigenvectors |
| `data/stats_seed/` | Cached category PC1 scores, retention map and extra city indices |
| `data/sewage/` | Derived sewage city indices and environmental feature tables |
| `data/pathogen_species22.csv` | Actual 22-species subset for the 2,049 included surface samples |
| `data/provenance.json` | Included input origins and transformations |
| `scripts/` | Five portable statistical steps |
| `reference_results/` | Unmodified historical output tables used for comparisons |
| `validation/` | Recorded test environment and validation receipt |
| `SHA256SUMS.txt` | Integrity checksums for distributed files |

The five steps implement: (1) city factor/category screens with 10 spatial axes,
5-axis sensitivity and factor-level sample re-testing; (2) category-level re-testing;
(3) exact block-order LMG and composite checks; (4) sewage follow-up and
185-pair concordance; (5) dominant-species sensitivity.
The spatial-axis routine is preserved from the source scripts: positive eigenvectors
of double-centered squared great-circle distances, **without a truncation threshold**.
The historical function name `dbmem` should not be read as a new method validation.
The historical `replicated` CSV field means nominal direction-prespecified one-sided
p < 0.05 under the original operational definition, not mechanistic validation.

S14 retains the exact historical 22-species list and raw relative-abundance basis;
it does not reclose the remaining taxa after mock-taxon removal. The historical list
includes *Vibrio natriegens*; inclusion in this operational index is not a claim that
every taxon is a demonstrated human pathogen. This package does not re-curate that list.
The documented Zurich station-latitude correction (47.0 degrees north) is retained.
Deterministic category sorting was added for LMG iteration; it does not change the
method, and comparisons ignore the text ordering of category names.

## Not covered

No claim is made that this archive regenerates every supplementary table, upstream
quality-control decision, raw taxonomic/AMR/pathway profile, environmental extraction,
or embedded manuscript figure. It does not provide independent replication.
No raw reads, large third-party databases, credentials, logs or private paths are
required by the distributed scripts.

## Provenance, permissions and citation

See `PROVENANCE.md`. Original data remain subject to their source terms.
No blanket open-source license is assigned here: the corresponding author must
confirm redistribution rights and choose a license for author-owned code before
public release. Repository creation or this ZIP alone does not make the materials
public. Add the final public repository URL and archived release identifier to the
manuscript availability statement after uploading. Cite the manuscript and original
source studies as described in the manuscript reference list; no DOI is invented.

中文上传步骤见 `UPLOAD_TO_GITHUB_ZH.md`。

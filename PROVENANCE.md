# Provenance and reproducibility boundary

The CSV inputs are real saved derived tables from the study workspace, not simulated
or synthesized data. `data/provenance.json` records their relative source names.
Reference outputs are copied unchanged from the final statistical output set and
the historical S14 dominance-sensitivity calculation. Scripts were adapted from
the corresponding historical scripts by replacing absolute paths with package-relative
locations; LMG category order was sorted for deterministic output. The extracted
S14 script excludes all historical Word/XLSX editing code.

Discovery source: MetaSUB public surface metagenomes and the study's saved derived
functional indices. Sewage source: Global Sewage Surveillance Initiative public
PanRes resistance-gene and mOTUs genus profiles, subsequently summarized as saved city
indices. Source-study citations and accession details remain in the manuscript;
the data and code deposited by the source studies should also be credited.
This package starts after those processing steps and does not independently verify
their reconstruction.

Environmental data families represented in the derived tables include atmospheric
reanalysis, satellite particulate matter, soil metals, river pharmaceutical
concentrations, pesticides, socioeconomic indicators, built environment, urban
heat islands, climate, and soil properties. The manuscript Methods and original
source references specify these resources. The included derived tables are not
substitutes for the original source products or their licenses.

The corresponding author should confirm that these derived tables can be redistributed
under the applicable original terms before making the repository public.
If some files cannot be redistributed, a rights-compliant replacement retrieval
procedure will be needed; the current self-contained run depends on these inputs.

Validation compares rerun statistics with saved outputs; this is computational
consistency on shared data, **not independent scientific validation**. Statistical
test families, covariates, transformations and historical construct mismatches are
retained rather than silently revised. See README for interpretation limits.

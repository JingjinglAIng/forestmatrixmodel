# forestmatrixmodel
Global forest matrix demographic model (MATRIX)
MATRIX Simulation (Unified) — README

This README accompanies matrix_simulation.R, a single, public-ready R script that runs MATRIX forest stand simulations in two modes:

NA_FT — North America workflow using Forest Type (FT) models and FT biomass vectors.

GLOBAL_GEC — Global (non-NA) workflow using continent-keyed models and Ecozone (GEC) biomass vectors.

The script discovers missing output chunks, supports SLURM arrays or CLI index, records annual AGB (total, upgrowth, recruitment, mortality), and optionally uses time-varying climate.

Quick Start
# GLOBAL mode (GEC biomass, models keyed by CONT):
Rscript matrix_simulation.R \
  --mode=GLOBAL_GEC \
  --input=/path/Grid3km_cov_Dvec_OneHot_500Chunk_0405.csv \
  --models=/path/models \
  --biomass=/path/AGB_by_GEZ_compiled_0228.csv \
  --out_dir=/path/25y_outputs \
  --clim=CC --years=25 --continent=Africa \
  --out_prefix=glob3km_onehot_25y

# North America mode (FT biomass, FT-specific model files):
Rscript matrix_simulation.R \
  --mode=NA_FT \
  --input=/path/Grid3km_cov_Dvec_OneHot_500Chunk_0411.csv \
  --models=/path/models \
  --biomass=/path/FT_biomass_kg_0227_Mitra.csv \
  --out_dir=/path/25y_outputs \
  --clim=CC --years=25 --continent=NAmerica \
  --out_prefix=glob3km_nam_25y


The script auto-selects the missing ChunkIDs (based on existing files in --out_dir). Provide the target index via SLURM_ARRAY_TASK_ID or the first CLI positional arg (1-based).

Installation & Requirements

R ≥ 4.1

R packages:

randomForest

stringr

Install packages:

install.packages(c("randomForest","stringr"))

Inputs
1) Mapping Grid CSV (--input)

A chunked grid containing stand structure, site covariates, and identifiers.

Required columns (by mode):

Common:

ID, ChunkID, CONTINENT (e.g., NAmerica, Africa), optional LAT, LON

DBH distribution: DBH1…DBH13 (trees/ha by DBH class)

Climate covariates (constant): C1…C21

Optional scenario columns: C1S1…C21S1 (RCP45), C1S2…C21S2 (RCP85)

NA_FT additionally needs: FT (forest type code)

GLOBAL_GEC additionally needs: GEC (ecozone code)

The script sets dY = --years internally.

2) Models directory (--models)

NA_FT expects exact FT-specific files:

RF.UG.0503.ft.<FT>.RData

RF.RC.0426.ft.<FT>.RData

RF.MT.0425.ft.<FT>.RData

GLOBAL_GEC expects continent-keyed files (first match is used):

RF.UG.*<CONT>.RData

RF.RC.*<CONT>.RData

RF.MT.*<CONT>.RData

Each file defines an object RF (a randomForest model).

3) Biomass CSV (--biomass)

NA_FT: must contain columns named FT<FT> with 13 values each (kg/tree per DBH class).

GLOBAL_GEC: must contain columns named by GEC (e.g., Boreal, TropMoist), each with 13 values.

The script selects the correct column per row and converts to Mg/ha (/ 1000).

4) Bioclimate CSV (--bioclim, optional)

A wide table with annual climate covariates indexed by ID. If omitted, the script uses constant climate from the covariates in --input.

Outputs

For each processed ChunkID, the script writes:

<out_dir>/<out_prefix>_<ChunkID>.csv


Schema (superset; some columns may be NA if not in input):

Core: ID, LAT, LON, GEC, FT, CONT, Year, B, N, S1, S2

Initial DBH distribution: TPH1_1…TPH1_13

Final DBH distribution: TPH2_1…TPH2_13

Annual sequence (for dY = --years):

Y1, AGB1, UP1, RC1, MT1, Y2, AGB2, UP2, RC2, MT2, ...

AGB = total aboveground biomass (Mg/ha)

UP, RC, MT = annual components (Mg/ha): upgrowth, recruitment, mortality

Year is the number of simulated years (e.g., 25). Internal annual timestamps use 1999 + m.

Command-Line Flags
Flag	Required	Default	Description
--mode	✳︎	GLOBAL_GEC	NA_FT or GLOBAL_GEC
--input	✳︎	—	Input grid CSV
--models	✳︎	—	Directory with RF models
--biomass	✳︎	—	Biomass conversion CSV
--out_dir	✳︎	—	Output directory
--out_prefix		glob3km_onehot_25y (GLOBAL) / glob3km_nam_25y (NA)	Output filename prefix
--clim		CC	CC, RCP45, RCP85
--years		25	Number of years to simulate (dY)
--continent		(none)	Filter input to a single CONTINENT (recommended for mode separation)
--bioclim		(none)	Bioclimate annual CSV; if omitted, constant climate used
--quiet		false	Reduce console messages

Indexing:

SLURM_ARRAY_TASK_ID or first positional arg (1-based) selects the nth missing chunk to compute.

SLURM Examples
# GLOBAL mode
sbatch --array=1-100 --wrap \
'Rscript matrix_simulation.R \
  --mode=GLOBAL_GEC \
  --input=/scratch/.../Grid3km_cov_Dvec_OneHot_500Chunk_0405.csv \
  --models=/scratch/.../models \
  --biomass=/scratch/.../AGB_by_GEZ_compiled_0228.csv \
  --out_dir=/scratch/.../25y_0820 \
  --clim=CC --years=25 --continent=Africa \
  --out_prefix=glob3km_onehot_25y --quiet=true'

# NA mode
sbatch --array=1-200 --wrap \
'Rscript matrix_simulation.R \
  --mode=NA_FT \
  --input=/scratch/.../Grid3km_cov_Dvec_OneHot_500Chunk_0411.csv \
  --models=/scratch/.../models \
  --biomass=/scratch/.../FT_biomass_kg_0227_Mitra.csv \
  --out_dir=/scratch/.../25y_0820 \
  --clim=CC --years=25 --continent=NAmerica \
  --out_prefix=glob3km_nam_25y --quiet=true'

How It Works (Brief)

Chunk discovery: Finds all ChunkID in --input, compares to files in --out_dir matching --out_prefix_*.csv, and identifies missing ones.

Model selection:

NA_FT: loads FT-specific RF.UG.0503.ft.<FT>.RData, RF.RC.0426.ft.<FT>.RData, RF.MT.0425.ft.<FT>.RData.

GLOBAL_GEC: picks the first file matching RF.UG.*<CONT>.RData (and likewise for RC/MT).

Annual loop: updates DBH class vector with upgrowth, stasis, mortality, and recruitment; computes AGB and components using the appropriate biomass vector (FT or GEC).

Output: writes one CSV per chunk; atomic write via *.tmp_$$ then rename, with a .lock_<prefix>_<ChunkID> to avoid duplicate work.

Data Conventions & Units

DBH classes: 13 bins; internal mean DBH vector is fixed in the script.

AGB and its components (UP, RC, MT) are in Mg/ha.

Diversity:

S1 = Shannon index

S2 = Simpson index

Biomass vectors are kg per tree per DBH class (converted to Mg/ha by /1000).

Troubleshooting

“Biomass vector not found”:
Ensure your biomass CSV has columns named:

NA_FT: FT<FT> (e.g., FTPine, FTMixedHardwood)

GLOBAL_GEC: names exactly matching GEC values in your input

“Cannot find model files for CONT” (GLOBAL_GEC):
Check that model filenames contain the continent string from CONTINENT (e.g., Africa, SAmerica, Europe).

All chunks processed but job still running:
The script exits early if nothing is missing; ensure --out_prefix matches what’s already written.

Parallel collisions:
Lock files .lock_<prefix>_<ChunkID> prevent duplicate work. If a crashed job leaves a stale lock, remove it only after confirming no running job uses that chunk.

Climate scenarios:
If --bioclim is omitted, the script uses constant climate:

--clim=CC: C1..C21

--clim=RCP45: C1S1..C21S1

--clim=RCP85: C1S2..C21S2

Reproducibility Tips

Record:

Git commit (if you version your inputs/models)

R version and package versions

Command line used (all flags)

Pin model filenames and biomass tables in a manifest.

Citation

If this work contributes to a publication, please cite the MATRIX modeling framework and related datasets/models from your Methods, and acknowledge:

Liang, J. et al. (under review) The hidden demography of the 21st century global forest carbon sink.

License
MIT

Maintainer

Jingjing Liang — albeca.liang@gmail.com

Contributions, issues, and improvements are welcome.

ChatGPT can make mistakes. OpenAI doesn't use FACAI-Purdue workspace data to train its models.

MATRIX Simulation (Unified) — README

This README accompanies matrix_simulation.R, a public-ready R script that runs MATRIX forest stand simulations in two modes:

NA_FT — North America workflow using Forest Type (FT) models and FT biomass vectors.

GLOBAL_GEC — Global (non-NA) workflow using continent-keyed models and Ecozone (GEC) biomass vectors.

The script discovers missing output chunks, supports SLURM arrays or CLI index, records annual AGB (total, upgrowth, recruitment, mortality), and optionally uses time-varying climate.

📦 All trained random forest models are provided in the models.zip file.
Unzip this file and point --models to the extracted directory.

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

install.packages(c("randomForest", "stringr"))

Inputs

Mapping Grid CSV (--input)
A chunked grid containing stand structure, site covariates, and identifiers.

Required columns (by mode):

Common:
ID, ChunkID, CONTINENT, optional LAT, LON
DBH1…DBH13 (trees/ha per DBH class)
Climate covariates: C1…C21
Optional: C1S1…C21S1 (RCP45), C1S2…C21S2 (RCP85)

NA_FT: requires FT (forest type code)

GLOBAL_GEC: requires GEC (ecozone code)

The script sets dY = --years internally.

Models directory (--models)
Unzip the provided models.zip archive before running.

NA_FT expects FT-specific files:

RF.UG.0503.ft.<FT>.RData

RF.RC.0426.ft.<FT>.RData

RF.MT.0425.ft.<FT>.RData

GLOBAL_GEC expects continent-keyed files (first match is used):

RF.UG.*<CONT>.RData

RF.RC.*<CONT>.RData

RF.MT.*<CONT>.RData

Each file defines an object RF (a randomForest model).

Biomass CSV (--biomass)

NA_FT: must contain columns named FT<FT> (13 values each, kg/tree by DBH class).

GLOBAL_GEC: must contain columns named by GEC (e.g., Boreal, TropMoist).

The script selects the correct column per row and converts to Mg/ha (/1000).

Bioclimate CSV (--bioclim, optional)
Annual climate covariates indexed by ID. If omitted, the script uses constant climate from --input.

Outputs

For each processed ChunkID, the script writes:

<out_dir>/<out_prefix>_<ChunkID>.csv


Schema (superset; some columns may be NA if not in input):

Core: ID, LAT, LON, GEC, FT, CONT, Year, B, N, S1, S2

Initial DBH distribution: TPH1_1…TPH1_13

Final DBH distribution: TPH2_1…TPH2_13

Annual sequence (for dY years):
Y1, AGB1, UP1, RC1, MT1, Y2, AGB2, ...

Units:

AGB, UP, RC, MT are in Mg/ha.

S1 = Shannon index, S2 = Simpson index.

Command-Line Flags
Flag	Required	Default	Description
--mode	✓	GLOBAL_GEC	Mode: NA_FT or GLOBAL_GEC
--input	✓	—	Input grid CSV
--models	✓	—	Directory with RF models (from models.zip)
--biomass	✓	—	Biomass conversion CSV
--out_dir	✓	—	Output directory
--out_prefix		Mode-dependent	Output filename prefix
--clim		CC	Climate scenario: CC, RCP45, RCP85
--years		25	Number of simulation years
--continent		—	Restrict to one CONTINENT
--bioclim		—	Bioclimate annual CSV (optional)
--quiet		false	Suppress console messages
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

Citation

If this work contributes to a publication, please cite the MATRIX modeling framework and related datasets/models from your Methods, and acknowledge:

Liang, J. et al. (under review). The hidden demography of the 21st century global forest carbon sink.

License

MIT License

Maintainer

Jingjing Liang — albeca.liang@gmail.com

Contributions, issues, and improvements are welcome.

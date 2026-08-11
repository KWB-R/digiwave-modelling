# ============================================================================
# R/config.R  —  single source of truth for ALL file paths + shared helpers
# ----------------------------------------------------------------------------
# Source this at the top of every script:
#     source("R/config.R")          # when the project root is the working dir
# (Opening DigiWave.Rproj in RStudio sets the working dir to the project root.)
#
# Every external input below resolves to a file under data/raw/. Drop the
# provided files there (same names) and every script reads them locally.
# ============================================================================

# ---- resolve project root (works in RStudio and via Rscript) ---------------
PROJECT_ROOT <- tryCatch(
  if (requireNamespace("here", quietly = TRUE)) here::here() else getwd(),
  error = function(e) getwd()
)

# ---- folder tree -----------------------------------------------------------
DATA_DIR    <- file.path(PROJECT_ROOT, "data")
RAW_DIR     <- file.path(DATA_DIR, "raw")       # <-- put provided external files here
DERIVED_DIR <- file.path(DATA_DIR, "derived")   # intermediate objects (.rds)
MODELS_DIR  <- file.path(PROJECT_ROOT, "models")
PLOTS_DIR   <- file.path(PROJECT_ROOT, "plots")
DOCS_DIR    <- file.path(PROJECT_ROOT, "docs")

for (d in c(RAW_DIR, DERIVED_DIR, MODELS_DIR, PLOTS_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- path helpers ----------------------------------------------------------
data_path    <- function(...) file.path(DATA_DIR,    ...)
raw          <- function(...) file.path(RAW_DIR,     ...)
derived_path <- function(...) file.path(DERIVED_DIR, ...)
model_path   <- function(...) file.path(MODELS_DIR,  ...)
plot_path    <- function(...) file.path(PLOTS_DIR,   ...)

# ---- shared-object persistence (makes each script standalone) --------------
# Producers call save_derived(obj, "name"); consumers call load_derived("name").
save_derived <- function(obj, name) {
  saveRDS(obj, derived_path(paste0(name, ".rds")))
  invisible(obj)
}
load_derived <- function(name) {
  f <- derived_path(paste0(name, ".rds"))
  if (!file.exists(f))
    stop("Missing intermediate object '", name, "' (", f, ").\n",
         "Run its producer script first — see the dependency table in README.md.",
         call. = FALSE)
  readRDS(f)
}
# Inject a saved list of objects into the calling environment by name.
attach_derived <- function(name, envir = parent.frame()) {
  list2env(load_derived(name), envir = envir)
}

# ============================================================================
# EXTERNAL INPUTS  (provide these under data/raw/)
# ============================================================================

## Master cleaned model input
MODEL_INPUT_RDATA <- raw("ModelInput_CLEAN_final.Rdata")

## aquaBio online E. coli analyser -------------------------------------------
AQUABIO_B403_20240417 <- raw("aquaBio_raw-data/20240417/historical_b403.csv")   # Aquabio.R
AQUABIO_B403_20240527 <- raw("aquaBio_raw-data/Datenformat_1/20240527/historical_b403.csv") # DataAnalysis.R / training.R
AQUABIO_ECOLI_2024_08 <- raw("aquaBio_raw-data/Datenformat_2/20240820/Historical_ecoli.csv") # DataAnalysis.R
AQUABIO_ECOLI_2024_06 <- raw("aquaBio_raw-data/Datenformat_2/20240605/Historical_ecoli.csv") # training.R
SBO_MIBI_NACH_GAK     <- raw("SBO_MiBi_Nach_GAK.csv")                            # Aquabio.R
# aquaBio post-UV history — already present locally in data/ (validation set)
AQUABIO_ECOLI_HIST    <- data_path("Historical_ecoli.csv")

## Plant (Klaeranlage) 15-min data -------------------------------------------
KA_15MIN_2024_08_19   <- raw("KA_15min_reshaped/20240819_Klaeranlage.csv")       # DataAnalysis.R
KA_15MIN_DATEN        <- raw("KA_15min_reshaped/Daten_Klaeranlage.csv")          # training.R

## Bench-scale UV dosing (dose-response) -------------------------------------
CBD_ECOLI_CLEAN2      <- raw("ecoli_clean2.csv")                                 # CBD.R
# UV dose vs LRV table used to fit stanDR
UV_DOSE_RESPONSE_CSV  <- raw("plot-data.csv")                                    # Laboratory.R

## Laboratory microbiology ----------------------------------------------------
LAB_MESSUNGEN_CSV     <- raw("20241111_DigiWave_Labormessungen.csv")             # Laboratory.R
ZULAUF_KA_CSV         <- raw("ZulaufKA.csv")                                     # preValidation.R
ABLAUF_UV_CSV         <- raw("AblaufUV.csv")                                     # preValidation.R

## Eurofins lab-report PDFs (folder of dated subfolders) ---------------------
EUROFINS_DIR          <- raw("Eurofins")                                         # read_pruefbericht.R

## Production prediction log (SQLite) + exported predictions ------------------
SINK_DB               <- raw("sink.prod.db")                                     # Database.R
PREDICTIONS_CSV       <- raw("predictions.csv")                                  # Laboratory.R

## Plant CSVs already present locally in data/ (referenced by name) ----------
KLAERWERK_S2025_CSV   <- data_path("KlaerwerkS2025.csv")                         # Database.R
KLAERWERK_NOVD_CSV    <- data_path("klaerwerkNOVD.csv")                          # CleanNewKAData.R (in_file)
UV_ZULAUF_CSV         <- data_path("uv_zulauf.csv")                              # several

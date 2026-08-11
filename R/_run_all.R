# ============================================================================
# _run_all.R  (PUBLIC copy) — reproduce the validation + energy results.
#
# Stages 01 (clean raw plant CSV) and 02 (train models) require the confidential
# WWTP sensor data and are therefore NOT runnable here — they are included as
# documented code only (raw data available from the authors on request).
#
# This copy ships the fitted models with their sensor split-thresholds
# affine-obscured, and the sensor features affine-transformed with the SAME map,
# so stages 03–06 reproduce the published figures. E. coli, timestamps and the
# UV dose-response are real. See README.md ("Data anonymization").
#
# Run from the project root:  Rscript R/_run_all.R
# ============================================================================

source("R/config.R")

pipeline <- c(
  "R/03_dose_response.R",     # -> models/stanDR.rds (from real dose-response data)
  "R/04_dose_schedule.R",     # -> data/derived/dose_schedule.rds
  "R/05_predict_validate.R",  # validation figures (uses shipped sanitized model + features)
  "R/06_energy_savings.R"     # dose recommendation + energy-saving figures
)

message("NOTE: stages 01_clean_plant_data.R and 02_train_models.R require the ",
        "confidential raw sensor data and are not run here (code-only).")

for (s in pipeline) {
  message("\n=== running ", s, " ===")
  ok <- tryCatch({ source(s, local = new.env()); TRUE },
                 error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) stop("Pipeline stopped at ", s, call. = FALSE)
}
message("\nPublic pipeline complete (stages 03–06).")

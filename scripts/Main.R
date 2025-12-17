# ======================================================
#                        MAIN RUNNER
# ======================================================
#   Executes both analytical pipelines:
#     1. GIS + Bayesian Analysis
#     2. Full Statistical Analysis
#
#   Author: Iñaki Intxaurbe Alberdi 
#   Department of Graphic Design and Engineering Projects
#   (Universidad del País Vasco/Euskal Herriko Unibertsitatea)
#   PACEA UMR 5199
#   (Université du Bordeaux)
#   Date: 2025-10-07
#   Copyright (C) 2025  Iñaki Intxaurbe
# ======================================================

# ---- 1. INITIAL SETUP ----
base_dir     <- dirname(getwd())  
data_dir     <- file.path(base_dir, "data")
scripts_dir  <- file.path(base_dir, "scripts")
outputs_dir  <- file.path(base_dir, "outputs")

# Create outputs directory if it doesn't exist
if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE)

# ---- 2. CHECK DATA FILE ----
data_file <- file.path(data_dir, "Table-Data-Base.xlsx")
if (!file.exists(data_file)) {
  stop("'Table-Data-Base.xlsx' not found in the /data/ folder. 
        Please make sure the file is correctly placed before running this script.")
} else {
  message(" Data file found at: ", data_file)
}

# ---- 3. FUNCTION TO RUN SCRIPTS ----
# Each script runs with working directory = /outputs/
run_script <- function(script_name, data_path, outputs_dir) {
  script_path <- file.path(scripts_dir, script_name)
  message("\n Running ", script_name, " ...")
  
  tryCatch({
    # create a local environment for isolated execution
    local_env <- new.env()
    local_env$data_path <- data_path
    
    # temporarily change working directory to outputs/
    old_wd <- getwd()
    setwd(outputs_dir)
    
    source(script_path, local = local_env)
    
    # restore working directory
    setwd(old_wd)
    
    message("Completed successfully: ", script_name)
  }, error = function(e) {
    message("Error in ", script_name, ": ", e$message)
  })
}

# ---- 4. EXECUTE SCRIPTS IN ORDER ----
run_script("00_GIS_and_Datings.R", data_file, outputs_dir)
run_script("01_Full_Statistics.R", data_file, outputs_dir)

# ---- 5. FINAL MESSAGE ----
message("\n ALL ANALYSES COMPLETED SUCCESSFULLY")
message("Results have been saved in the /outputs/ folder:")
message(" - 'Datings_results.xlsx' + GIS and Bayesian plots")
message(" - 'evidences_Tests.xlsx' + statistical plots")

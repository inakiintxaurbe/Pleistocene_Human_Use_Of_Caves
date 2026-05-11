#   MAIN RUNNER : Pleistocene human use of caves
# 
#   Author: Iñaki Intxaurbe Alberdi 
#   Copyright (C) 2025  Iñaki Intxaurbe
#
#   SPDX-License-Identifier: AGPL-3.0 (citation mandatory)
# ======================================================

base_dir     <- dirname(getwd())  
data_dir     <- file.path(base_dir, "data")
scripts_dir  <- file.path(base_dir, "scripts")
outputs_dir  <- file.path(base_dir, "outputs")

dir.create(outputs_dir, recursive = TRUE)

data_file <- file.path(data_dir, "Table-Data-Base.xlsx")

run_script <- function(script_name, data_path, outputs_dir) {
  script_path <- file.path(scripts_dir, script_name)
  
  local_env <- new.env()
  local_env$data_path <- data_path
    
  old_wd <- getwd()
  setwd(outputs_dir)
    
  source(script_path, local = local_env)
  setwd(old_wd)
}

run_script("00_GIS_and_Datings.R", data_file, outputs_dir)
run_script("01_Full_Statistics.R", data_file, outputs_dir)
run_script("02_Full_Statistics_Bayessian_Only.R", data_file, outputs_dir)

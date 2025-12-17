```
# Pleistocene_Human_Use_Of_Caves

Created by: Iñaki Intxaurbe Alberdi
Department of Graphic Design and Engineering Projects
(Universidad del País Vasco/Euskal Herriko Unibertsitatea)
PACEA UMR 5199
(Université du Bordeaux)
Date: 2025-10-07

Copyright (C) 2025  Iñaki Intxaurbe

This repository contains R scripts and data workflows developed to analyze the archaeological evidence of human use of caves during the Pleistocene. The project combines quantitative statistical methods, correspondence and multivariate analyses, and data visualization to explore relationships between cultural evidences, cave depth, and accessibility.

## Overview

The main analysis script (`01_Full_Statistics.R`) processes a standardized archaeological database (`Table-Data-Base.xlsx`) and produces:
- Statistical summaries (Chi-square, Kruskal–Wallis, residual analyses)
- Multivariate analyses (CA, MCA)
- Visualization outputs (heatmaps, boxplots, violin plots)
- An Excel report summarizing all test results (`evidences_Tests.xlsx`)

## Folder Structure

Pleistocene_Human_Use_Of_Caves/
│
├── 01_Full_Statistics.R # Main R analysis script
├── Table-Data-Base.xlsx # Input archaeological dataset
├── Plots/ # Automatically generated figures
│ ├── Heatmap_* # Heatmaps for evidences, depth, etc.
│ ├── Boxplot_* # Boxplots by region and phase
│ └── Violin_* # Violin plots for depth distributions
└── evidences_Tests.xlsx # Statistical output summary

## Requirements

Before running the analysis, ensure that R (≥ 4.1) and the following packages are installed:

```r
install.packages(c(
  "readxl", "dplyr", "FactoMineR", "factoextra", "pheatmap",
  "openxlsx", "reshape2", "tibble", "tidyr", "stringr",
  "ggplot2", "DescTools", "RColorBrewer"
))
```
### **Running the Analysis**

Place the Excel file _Table-Data-Base.xlsx_ in the project root directory.

Open _01_Full_Statistics.R_ in RStudio (or your preferred R environment).

Run the script. It will:

Perform all statistical computations

Generate visual outputs in the Plots/ folder

Save results in evidences_Tests.xlsx

At the end of execution, a message will confirm completion:

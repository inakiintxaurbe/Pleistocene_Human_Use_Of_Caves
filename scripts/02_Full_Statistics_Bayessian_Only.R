#   FULL STATISTICAL ANALYSIS (BAYESIAN ONLY) 
# 
#   Author: Iñaki Intxaurbe Alberdi 
#   Copyright (C) 2025  Iñaki Intxaurbe
#
#   SPDX-License-Identifier: AGPL-3.0 (citation mandatory)

# Install packages and organise framework ----------------------------------------------------------------------------------------------------------------------------------------------

pkgs <- c(
  "readxl","dplyr","FactoMineR","factoextra",
  "pheatmap","openxlsx","reshape2","tibble","tidyr",
  "stringr","ggplot2","DescTools","RColorBrewer"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(readxl)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(openxlsx)
library(reshape2)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)
library(RColorBrewer)

# Helpers / Laguntzailiak---------------------------------

roman_levels <- c("I","II","III","IV","V","VI","VII","VIII","IX","X")

roman2int <- function(x) suppressWarnings(as.numeric(as.roman(trimws(toupper(x)))))

theme_pub <- theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ------------------ DATA INPUT -------------------------

file_path <- if (exists("data_path")) data_path else file.path(getwd(), "Table-Data-Base.xlsx")

# 1) CAVES (approved only)
df_caves <- read_excel(file_path, sheet = "Caves") %>%
  filter(tolower(Approved) == "yes")

# 2) BAYESIAN PHASES (from Monte-Carlo results)
df_bayes <- read_excel("Datings_results.xlsx", sheet = "Datings_combined") %>%
  distinct(Cave, Phase)

# 3) KEEP ONLY BAYESIAN CAVE–PHASE PAIRS
df <- df_caves %>%
  inner_join(df_bayes, by = c("Cave","Phase"))


# ------------------ EVIDENCES --------------------------

evidences <- c(
  "Rock Art","Portable Art","Structures / Speleofacts",
  "Lithic Industry","Modified bones","Ochre remains",
  "Human Remains","Footprints","Fire Remains"
)

df_bin <- df %>%
  mutate(across(all_of(evidences), ~{
    x <- toupper(trimws(as.character(.)))
    as.integer(ifelse(!is.na(x) & x == "X", 1, 0))
  })) %>%
  mutate(
    Phase_clean = factor(Phase, levels = roman_levels, ordered = TRUE)
  )


# ------------------ DEPTH & DIFFICULTY -----------------

# Depth
df_bin <- df_bin %>%
  mutate(Depth_m = suppressWarnings(
    as.numeric(gsub(",", ".", as.character(`Reached aprox. depth (m)`)))
  ))

# Difficulty
df_bin <- df_bin %>%
  mutate(
    Difficulty_num = case_when(
      `Level of difficulty of access` == "Very Low"  ~ 1,
      `Level of difficulty of access` == "Low"       ~ 2,
      `Level of difficulty of access` == "Medium"    ~ 3,
      `Level of difficulty of access` == "Hard"      ~ 4,
      `Level of difficulty of access` == "Very Hard" ~ 5,
      TRUE ~ NA_real_
    ),
    Difficulty_fac = factor(
      Difficulty_num,
      levels = 1:5,
      labels = c("Very Low","Low","Medium","Hard","Very Hard")
    )
  )


# ------------------ OUTPUT FOLDER ----------------------

plot_dir <- file.path(getwd(), "Plots_Bayesian_Only")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# ======================================================
#                     HEATMAPS
# ======================================================

## Phase × Evidence (Bayesian only)
mat_phase <- df_bin %>%
  group_by(Phase_clean) %>%
  summarise(across(all_of(evidences), ~ mean(.==1, na.rm=TRUE)), .groups="drop") %>%
  as.data.frame()

rownames(mat_phase) <- mat_phase$Phase_clean
mat_phase <- mat_phase[, -1, drop=FALSE]

png(file.path(plot_dir,"Heatmap_Evidences_Phase_Bayesian.png"),
    width=2400, height=1800, res=300)

pheatmap(as.matrix(mat_phase),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 10,
         angle_col = 45,
         main = "Evidences by Phase (Bayesian only)")

dev.off()


## Region × Evidence (Bayesian only)
mat_region <- df_bin %>%
  group_by(Region) %>%
  summarise(across(all_of(evidences), ~ mean(.==1, na.rm=TRUE)), .groups="drop") %>%
  column_to_rownames("Region")

png(file.path(plot_dir,"Heatmap_Evidences_Region_Bayesian.png"),
    width=2400, height=1800, res=300)

pheatmap(as.matrix(mat_region),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 10,
         angle_col = 45,
         main = "Evidences by Region (Bayesian only)")

dev.off()


# ======================================================
#                     BOXPLOTS
# ======================================================

## Depth by Phase
if (!all(is.na(df_bin$Depth_m))) {
  g1 <- ggplot(df_bin, aes(x = Phase_clean, y = Depth_m, fill = Phase_clean)) +
    geom_boxplot() +
    labs(title="Depth by Phase (Bayesian only)",
         x="Phase", y="Depth (m)") +
    theme_pub
  
  ggsave(file.path(plot_dir,"Boxplot_Depth_Phase_Bayesian.png"),
         g1, width=9, height=6, dpi=300)
}

## Difficulty by Phase
if (!all(is.na(df_bin$Difficulty_num))) {
  g2 <- ggplot(df_bin, aes(x = Phase_clean, y = Difficulty_num, fill = Phase_clean)) +
    geom_boxplot() +
    scale_y_continuous(breaks=1:5,
                       labels=c("Very Low","Low","Medium","Hard","Very Hard")) +
    labs(title="Difficulty by Phase (Bayesian only)",
         x="Phase", y="Difficulty") +
    theme_pub
  
  ggsave(file.path(plot_dir,"Boxplot_Difficulty_Phase_Bayesian.png"),
         g2, width=9, height=6, dpi=300)
}


# ======================================================
#               CORRESPONDENCE ANALYSIS
# ======================================================

## CA – Phase × Evidence
tab_ca_phase <- df_bin %>%
  group_by(Phase_clean) %>%
  summarise(across(all_of(evidences), ~ sum(.==1, na.rm=TRUE)), .groups="drop") %>%
  as.data.frame()

rownames(tab_ca_phase) <- tab_ca_phase$Phase_clean
tab_ca_phase <- tab_ca_phase[, -1, drop=FALSE]
tab_ca_phase <- tab_ca_phase[rowSums(tab_ca_phase)>0, colSums(tab_ca_phase)>0]

if (nrow(tab_ca_phase) > 1 && ncol(tab_ca_phase) > 1) {
  ca_ph <- CA(tab_ca_phase, graph=FALSE)
  
  pdf(file.path(plot_dir,"CA_Phases_Evidences_Bayesian.pdf"), width=9, height=7)
  print(fviz_ca_biplot(ca_ph, repel=TRUE, ggtheme=theme_minimal()))
  dev.off()
}

### VIOLIN PLOTS / BIBOLIN FORMAKO GRAFIKUAK -------------

# Sakonera / Depth

if (!all(is.na(df_bin$Depth_m))) {
  g_violin_depth <- ggplot(
    df_bin %>% filter(!is.na(Depth_m)),
    aes(x = Phase_clean, y = Depth_m, fill = Phase_clean)
  ) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
    labs(
      title = "Depth by Phase (Bayesian only)",
      x = "Phase",
      y = "Reached aprox. depth (m)"
    ) +
    theme_pub
  
  ggsave(
    file.path(plot_dir, "Violin_Depth_Phase_Bayesian.png"),
    g_violin_depth,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# Zailtasuna / Difficulty ------------------------------

if (!all(is.na(df_bin$Difficulty_num))) {
  g_violin_diff <- ggplot(
    df_bin %>% filter(!is.na(Difficulty_num)),
    aes(x = Phase_clean, y = Difficulty_num, fill = Phase_clean)
  ) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
    scale_y_continuous(
      breaks = 1:5,
      labels = c("Very Low","Low","Medium","Hard","Very Hard")
    ) +
    labs(
      title = "Difficulty of access by Phase (Bayesian only)",
      x = "Phase",
      y = "Level of difficulty of access"
    ) +
    theme_pub
  
  ggsave(
    file.path(plot_dir, "Violin_Difficulty_Phase_Bayesian.png"),
    g_violin_diff,
    width = 9,
    height = 6,
    dpi = 300
  )
}


# ======================================================
#                    END
# ======================================================

message("Bayesian-only statistical analysis completed.")

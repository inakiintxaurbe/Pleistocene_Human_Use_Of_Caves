#   FULL STATISTICAL ANALYSIS
# 
#   Author: Iñaki Intxaurbe Alberdi 
#   Copyright (C) 2025  Iñaki Intxaurbe
#
#   SPDX-License-Identifier: AGPL-3.0 (citation mandatory)

# Install packages and organise framework ----------------------------------------------------------------------------------------------------------------------------------------------

pkgs <- c(
  "readxl","dplyr","FactoMineR","factoextra",
  "pheatmap","openxlsx","reshape2","tibble","tidyr",
  "stringr","ggplot2","DescTools","viridisLite"
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
library(viridisLite)

roman2int <- function(x) suppressWarnings(as.numeric(as.roman(trimws(toupper(x)))))

theme_pub <- theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

roman_levels <- c(
  "I",
  "II",
  "III",
  "IV",
  "V",
  "VI",
  "VII",
  "VIII",
  "IX",
  "X"
)

effect_pos <- function(extra="") paste0("Positive (more than expected", ifelse(extra=="","",paste0("; ",extra)), ")")
effect_neg <- function(extra="") paste0("Negative (less than expected", ifelse(extra=="","",paste0("; ",extra)), ")")
effect_none <- function() "No clear association"

project_dir <- normalizePath(
  file.path(
    getwd(), 
    ".."
  ), 
  winslash = "/"
)

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
dir.create(
  output_dir, 
  recursive = TRUE, 
  showWarnings = FALSE
)
plot_dir <- file.path(output_dir, "Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
save_plot <- function(
    filename, 
    plot, 
    w=10, 
    h=6, 
    dpi=300
  ) 
  {
  ggplot2::ggsave(
    filename = file.path(plot_dir, filename), 
    plot = plot, 
    width = w, 
    height = h, 
    dpi = dpi
  )
}

file_path <- file.path(data_dir, "Table-Data-Base.xlsx")
df <- read_excel(file_path, sheet = "Caves") %>% filter(tolower(Approved) == "yes")

evidences <- c(
  "Rock Art",
  "Portable Art",
  "Structures / Speleofacts",
  "Lithic Industry",
  "Modified bones",
  "Ochre remains",
  "Human Remains",
  "Footprints",
  "Fire Remains"
)

df_bin <- df %>%
  mutate(across(all_of(evidences), ~{
    x <- toupper(trimws(as.character(.)))
    as.integer(ifelse(!is.na(x) & x == "X", 1, 0))
  })) %>%
  mutate(
    Phase_clean = ifelse(grepl(
      "^(I|II|III|IV|V|VI|VII|VIII|IX|X)$", 
      trimws(Phase)), 
      trimws(Phase), 
      NA
    ),
    Phase_clean = factor(
      Phase_clean, 
      levels = roman_levels, 
      ordered = TRUE
    )
  )

region_levels <- sort(unique(df_bin$Region))

region_pal <- setNames(
  viridisLite::viridis(
    length(region_levels),
    option = "viridis"
  ),
  region_levels
)

phase_pal  <- setNames(
  viridisLite::viridis(
    length(roman_levels),
    option = "plasma"
  ),
  roman_levels
)


# Chi-square + residuals (evidences) ----------------------------------------------------------------------------------------------------------------------------------------------


results <- tibble(
  Evidence=character(), 
  Test=character(), 
  p_value=double(),
                  
  Significant=character(), 
  Max_category=character()
)
detailed_residuals <- list()

run_chi2_block <- function(
    var_bin_name, 
    group_var, 
    df_data
  ) {
  tab <- table(df_data[[group_var]], df_data[[var_bin_name]])
  if (!all(dim(tab) > 1)) return(NULL)
  test <- tryCatch(chisq.test(tab), error=function(e) fisher.test(tab))
  signif <- ifelse(test$p.value < 0.05, "Significant", "No Significant")
  cat_max <- NA
  if ("residuals" %in% names(test)) {
    res <- as.data.frame(as.table(test$residuals))
    colnames(res) <- c("Category","Value","Residual")
    res$Evidence <- var_bin_name; res$Test <- ifelse(group_var=="Phase_clean","Phase","Region")
    detailed_residuals[[length(detailed_residuals)+1]] <<- res
    cat_max <- res$Category[which.max(abs(res$Residual))]
  }
  tibble(Evidence=var_bin_name, Test=ifelse(group_var=="Phase_clean","Phase","Region"),
         p_value=test$p.value, Significant=signif, Max_category=cat_max)
}

for (v in evidences) {
  x1 <- run_chi2_block(
        v, 
        "Region", 
        df_bin
      );      if (!is.null(x1)) results <- bind_rows(results, x1)
  x2 <- run_chi2_block(
        v, 
        "Phase_clean", 
        df_bin
      ); if (!is.null(x2)) results <- bind_rows(results, x2)
}

detailed_residuals <- if (length(detailed_residuals)>0) 
  bind_rows(detailed_residuals) else tibble()

if (nrow(detailed_residuals) > 0) {
  detailed_residuals <- detailed_residuals %>%
    mutate(
      Value = suppressWarnings(as.integer(as.character(Value))),
      Its_presence = (Value == 1),
      Direction = case_when(
        Residual >=  2 ~ "Positive (more than expected)",
        Residual <= -2 ~ "Negative (less than expected)",
        TRUE ~ "No clear association"
      ),
      Es_Significant = ifelse(abs(Residual) >= 2, "Yes", "No")
    )
}

interpretations <- if (nrow(detailed_residuals) > 0) {
  cols <- c("Evidence", "Test", "Category", "Residual", "Direction")
  cols <- cols[cols %in% names(detailed_residuals)]
  detailed_residuals %>%
    filter(Its_presence, abs(Residual) >= 2) %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    arrange(dplyr::across(all_of(intersect(c("Test", "Evidence"), names(.)))), desc(Residual))
} else tibble()


# Clean depth & difficulty ----------------------------------------------------------------------------------------------------------------------------------------------

depth_col <- "Reached aprox. depth (m)"
stopifnot(depth_col %in% names(df_bin))
df_bin <- df_bin %>%
  mutate(Depth_m = suppressWarnings(
    as.numeric(
      gsub(
        ",", 
        ".", 
        as.character(.data[[depth_col]])
        )
      )
    )
  )

diff_col <- "Level of difficulty of access"
stopifnot(diff_col %in% names(df_bin))
normalize_difficulty <- function(x) {
  x0 <- trimws(as.character(x))  # cleans spaces at start/end
  
  case_when(
    x0 == "Very Low"  ~ 1,
    x0 == "Low"       ~ 2,
    x0 == "Medium"    ~ 3,
    x0 == "Hard"      ~ 4,
    x0 == "Very Hard" ~ 5,
    TRUE ~ NA_real_   # all other options are discarded
  )
}
df_bin <- df_bin %>%
  mutate(
    Difficulty_num = normalize_difficulty(.data[[diff_col]]),
    Difficulty_fac = factor(
                      Difficulty_num, 
                      levels = 1:5,
                      labels = c(
                          "Very Low",
                          "Low",
                          "Medium",
                          "Hard",
                          "Very Hard"
                        )
                      )
  )


# Tests: KW + CHI2 (Difficulty) + Residuals (Depth) ----------------------------------------------------------------------------------------------------------------------------------------------

# Kruskal-Wallis with interpretation
                   
kw_results_depth <- tibble(
  Variable=character(), Test=character(), Grouping=character(),
  p_value=double(), n_groups=integer(), Significant=character(),
  Max_category=character(), Direction=character()
)
kw_results_diff  <- kw_results_depth[0,]

                   
# Fixed function without tidy eval; works with strings
                   
add_kw <- function(data, value_col, group_col, var_name) {
  # value_col y group_col deben ser nombres de columnas (texto)
  if (!(group_col %in% names(data)) || !(value_col %in% names(data))) return(NULL)
  
  dd <- data %>%
    dplyr::select(all_of(c(group_col, value_col))) %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]]))
  
  if (nrow(dd) == 0 || length(unique(dd[[group_col]])) < 2) return(NULL)

  # Kruskal-Wallis test / KW test-a
  
  f <- as.formula(paste(value_col, "~", group_col))
  fm <- suppressWarnings(stats::kruskal.test(f, data = dd))

  # Median per group / Mediana taldeka
  
  med_by_grp <- tapply(dd[[value_col]], dd[[group_col]], median, na.rm = TRUE)
  max_cat <- names(which.max(med_by_grp))
  min_cat <- names(which.min(med_by_grp))
  
  # Interpretation / Interpretazioa
  
  direction <- if (fm$p.value < 0.05) {
    effect_pos(paste0("higher median in ", max_cat, " vs ", min_cat))
  } else effect_none()
  
  tibble(
    Variable = var_name,
    Test = "Kruskal-Wallis",
    Grouping = group_col,
    p_value = unname(fm$p.value),
    n_groups = length(med_by_grp),
    Significant = ifelse(fm$p.value < 0.05, "Significant", "No Significant"),
    Max_category = max_cat,
    Direction = direction
  )
}
                   
# Run Kruskal-Wallis tests / KW test-a burutu -------------------------

kw_results_depth <- bind_rows(
  kw_results_depth,
  add_kw(
    df_bin, 
    "Depth_m", 
    "Region", 
    "Depth_m"
  ),
  add_kw(
    df_bin, 
    "Depth_m", 
    "Phase_clean", 
    "Depth_m"
  )
)

kw_results_diff <- bind_rows(
  kw_results_diff,
  add_kw(
    df_bin, 
    "Difficulty_num", 
    "Region", 
    "Difficulty_num"
    ),
  add_kw(
    df_bin, 
    "Difficulty_num", 
    "Phase_clean", 
    "Difficulty_num"
    )
)

                   
# Chi² Difficulty (factor) with Fisher fallback + residuals ----------------------------------------------------------------------------------------------------------------------------------------------

assoc_results <- tibble(
  Variable=character(), Test=character(), Grouping=character(),
  p_value=double(), Assoc=character(),
  Significant=character(), Max_category=character(), Direction=character()
)

detailed_residuals_diff <- list()

run_diff_chi2 <- function(group_var) {
  # Filter valid data
  x <- df_bin %>% filter(!is.na(Difficulty_fac), !is.na(.data[[group_var]]))
  if (nrow(x) == 0 || length(unique(x[[group_var]])) < 2) return(NULL)
  
  # Force fixed levels
  diff_levels_fixed <- c("Very Low","Low","Medium","Hard","Very Hard")
  rows_levels <- sort(unique(na.omit(df_bin[[group_var]])))
  
  # Full table
  tab <- table(
    factor(x[[group_var]], levels = rows_levels),
    factor(as.character(x$Difficulty_fac), levels = diff_levels_fixed)
  )
  
# Remove empty rows/columns before testing 
  
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    tab_test <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab_test <- tab_test[, colSums(tab_test) > 0, drop = FALSE]
  } else {
    tab_test <- tab
  }
  if (!all(dim(tab_test) > 1)) return(NULL)
  
 # Chi² test with fallback to Fisher
  
  test <- tryCatch(suppressWarnings(chisq.test(tab_test)), error=function(e) fisher.test(tab_test))
  
  # Build complete residual table
  diff_levels_fixed <- c(
    "Very Low",
    "Low",
    "Medium",
    "Hard",
    "Very Hard"
  )
  rows_levels <- sort(unique(na.omit(df_bin[[group_var]])))
  full_grid <- expand.grid(Category = rows_levels, Value = diff_levels_fixed)
  
  if ("residuals" %in% names(test)) {
    res_tbl <- as.data.frame(as.table(test$residuals))
    names(res_tbl) <- c("Category","Value","Residual")
    
    full_grid <- full_grid %>%
      mutate(Category = as.character(Category), Value = as.character(Value))
    res_tbl <- res_tbl %>%
      mutate(Category = as.character(Category), Value = as.character(Value))
    
    res_tbl <- full_grid %>%
      left_join(
        res_tbl, 
        by = c("Category","Value")
      ) %>%
      mutate(
        Evidence = "Difficulty",
        Test = ifelse(
          group_var=="Phase_clean",
          "Phase",
          "Region"
        ),
        Direction = case_when(
          Residual >=  2  ~ "Positive (more than expected)",
          Residual <= -2  ~ "Negative (less presence than expected)",
          TRUE            ~ "No clear association"
        ),
        Es_Significant = ifelse(
          !is.na(Residual) & abs(Residual) >= 2, 
          "Yes", 
          "No"
        )
      )
  } else {
    full_grid <- full_grid %>%
      mutate(
        Category = as.character(Category), 
        Value = as.character(Value)
      )
    res_tbl <- full_grid %>%
      mutate(
        Residual = NA, 
        Evidence = "Difficulty",
        Test = ifelse(
          group_var=="Phase_clean",
          "Phase",
          "Region"
        ),
        Direction = "No data", Es_Significant = "No"
      )
  }
  
  detailed_residuals_diff[[length(detailed_residuals_diff)+1]] <<- res_tbl
                   
  # General direction
                   
  direction <- effect_none()
  if (exists("res_tbl") && nrow(res_tbl) > 0 && any(!is.na(res_tbl$Residual))) {
    idx <- which.max(abs(res_tbl$Residual))
    if (!is.na(res_tbl$Residual[idx])) {
      direction <- if (res_tbl$Residual[idx] >= 2) {
        effect_pos(paste0("peak in ", res_tbl$Category[idx], " for '", as.character(res_tbl$Value[idx]), "'"))
      } else if (res_tbl$Residual[idx] <= -2) {
        effect_neg(paste0("deficit in ", res_tbl$Category[idx], " for '", as.character(res_tbl$Value[idx]), "'"))
      } else effect_none()
    }
  }
                  
  # Cramer's V
                   
  cv <- tryCatch(DescTools::CramerV(tab_test, conf.level = NA), error=function(e) NA)
  
  tibble(
    Variable = "Difficulty_fac",
    Test = "Chi-square",
    Grouping = ifelse(
      group_var=="Phase_clean",
      "Phase",
      "Region"
    ),
    p_value = unname(test$p.value),
    Assoc = ifelse(
      is.na(cv), 
      NA, 
      sprintf(
        "Cramer's V = %.3f", 
        cv
      )
    ),
    Significant = ifelse(
      test$p.value < 0.05, 
      "Significant", 
      "No Significant"
    ),
    Max_category = rownames(tab_test)[which.max(rowSums(tab_test))],
    Direction = direction
  )
}

assoc_results <- bind_rows(
  run_diff_chi2("Region"), 
  run_diff_chi2("Phase_clean")
)
detailed_residuals_diff <- if (length(detailed_residuals_diff)>0) bind_rows(detailed_residuals_diff) else tibble()
  
# Residuals for DEPTH (z-score per group) ----------------------------------------------------------------------------------------------------------------------------------------------
                 
mk_depth_residuals <- function(group_var, group_label) {
  dat <- df_bin %>% filter(
    !is.na(Depth_m), 
    !is.na(.data[[group_var]])
  )
  if (nrow(dat) == 0) return(NULL)
  m_grp <- tapply(
    dat$Depth_m, 
    dat[[group_var]], 
    mean, 
    na.rm=TRUE
  )
  m_all <- mean(
    dat$Depth_m, 
    na.rm=TRUE
  )
  sd_all <- sd(
    dat$Depth_m, 
    na.rm=TRUE
  )
  if (is.na(sd_all) || sd_all == 0) return(NULL)
  z_grp <- (m_grp - m_all) / sd_all
  tibble(
    Category = names(z_grp),
    Value    = "Mean Depth",
    Residual = as.numeric(z_grp),
    Evidence = "Depth",
    Test     = group_label
  ) %>%
    mutate(
      Direction = case_when(
        Residual >=  2 ~ "Positive (more than expected)",
        Residual <= -2 ~ "Negative (less presence than expected)",
        TRUE ~ "No clear association"
      ),
      Es_Significant = ifelse(abs(Residual) >= 2, "Yes", "No")
    )
}
detailed_residuals_depth <- bind_rows(
  mk_depth_residuals(
    "Region", 
    "Region"
  ),
  mk_depth_residuals(
    "Phase_clean", 
    "Phase"
  )
)
                 
# Interpretations (only significant residuals) ----------------------------------------------------------------------------------------------------------------------------------------------
                 
interpretations_difficulty <- if (nrow(detailed_residuals_diff) > 0) {
  cols <- c(
    "Evidence", 
    "Test", 
    "Category", 
    "Residual", 
    "Direction"
  )
  cols <- cols[cols %in% names(detailed_residuals_diff)]
  detailed_residuals_diff %>%
    filter(Es_Significant == "Yes") %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    arrange(
      dplyr::across(
        all_of(
          intersect(
            c("Test"), 
            names(.)
          )
        )
      ), 
      desc(Residual)
    )
} else tibble()

interpretations_depth <- if (nrow(detailed_residuals_depth) > 0) {
  cols <- c(
    "Evidence", 
    "Test", 
    "Category",
    "Residual", 
    "Direction"
  )
  cols <- cols[cols %in% names(detailed_residuals_depth)]
  detailed_residuals_depth %>%
    filter(Es_Significant == "Yes") %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    arrange(
      dplyr::across(
        all_of(
          intersect(
            c("Test"), 
            names(.)
          )
        )
      ), 
      desc(Residual)
    )
} else tibble()

                 
# Create Excel ----------------------------------------------------------------------------------------------------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb,"Chi2_Results")
addWorksheet(wb,"detailed_residuals")
addWorksheet(wb,"Interpretations_presence")
addWorksheet(wb,"Difficulty_Chi2")
addWorksheet(wb,"Difficulty_KW")
addWorksheet(wb,"Residuals_Difficulty")
addWorksheet(wb,"Interpretations_Difficulty")
addWorksheet(wb,"Depth_KW")
addWorksheet(wb,"Residuals_Depth")
addWorksheet(wb,"Interpretations_Depth")
                 
safe_select <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(data.frame()) 
  df %>% dplyr::select(dplyr::all_of(cols))
}

writeData(
  wb, 
  "Chi2_Results", 
  results
)

writeData(
  wb, 
  "detailed_residuals",
  safe_select(
    detailed_residuals,
    c(
      "Category",
      "Value",
      "Residual",
      "Evidence",
      "Test",
      "Direction",
      "Es_Significant"
    )
  )
)

writeData(
  wb, 
  "Interpretations_presence", 
  interpretations
)
writeData(
  wb, 
  "Difficulty_Chi2", 
  assoc_results
)
writeData(
  wb, 
  "Difficulty_KW", 
  kw_results_diff
)

writeData(
  wb, 
  "Residuals_Difficulty",
  safe_select(
    detailed_residuals_diff,
    c(
      "Category",
      "Value",
      "Residual",
      "Evidence",
      "Test",
      "Direction",
      "Es_Significant"
    )
  )
)

writeData(
  wb, 
  "Interpretations_Difficulty", 
  interpretations_difficulty
)
writeData(
  wb, 
  "Depth_KW", 
  kw_results_depth
)

writeData(
  wb, 
  "Residuals_Depth",
  safe_select(
    detailed_residuals_depth,
    c(
      "Category",
      "Value",
      "Residual",
      "Evidence",
      "Test",
      "Direction",
      "Es_Significant"
    )
  )
)

writeData(
  wb, 
  "Interpretations_Depth", 
  interpretations_depth
)

saveWorkbook(
  wb,
  file.path(output_dir, "evidences_Tests.xlsx"),
  overwrite = TRUE
)

# Heatmaps (orignals) ----------------------------------------------------------------------------------------------------------------------------------------------
                 
## Proportions Region x Evidence
                 
mat_region <- df_bin %>%
  group_by(Region) %>%
  summarise(
    across(
      all_of(evidences), 
      ~ mean(
        .==1, 
        na.rm=TRUE
      )
    ), 
    .groups="drop"
  ) %>%
  column_to_rownames("Region")
if (nrow(mat_region) > 0) {
  png(
    file.path(
      plot_dir,
      "Heatmap_evidences_Region.png"
    ), 
    width=2400, 
    height=1800, 
    res=300
  )
  pheatmap(
    as.matrix(mat_region), 
    cluster_rows=TRUE, 
    cluster_cols=TRUE,
    show_rownames=TRUE, 
    show_colnames=TRUE, 
    fontsize_row=12, 
    fontsize_col=10,
    angle_col=45, 
    main="Proportion of evidences by region"
  )
  dev.off()
}
                 
## Proportions Phase x Evidence
                 
mat_phase <- df_bin %>%
  filter(!is.na(Phase_clean)) %>%
  group_by(Phase_clean) %>%
  summarise(
    across(
      all_of(evidences), 
      ~ mean(
        .==1, 
        na.rm=TRUE)
    ), 
    .groups="drop"
  ) %>%
  mutate(
    order_num = as.integer(Phase_clean)
  ) %>% arrange(order_num) %>%
  as.data.frame()
rownames(mat_phase) <- mat_phase$Phase_clean
mat_phase <- mat_phase %>% dplyr::select(
  -Phase_clean, 
  -order_num
)
if (nrow(mat_phase) > 0) {
  png(
    file.path(
      plot_dir,
      "Heatmap_evidences_Phase.png"
    ), 
    width=2400, 
    height=1800, 
    res=300
  )
  pheatmap(
    as.matrix(mat_phase), 
    cluster_rows=FALSE, 
    cluster_cols=TRUE,
    show_rownames=TRUE, 
    show_colnames=TRUE, 
    fontsize_row=12, 
    fontsize_col=10,
    angle_col=45, 
    main="Proportion of evidences by phase"
  )
  dev.off()
}


# Residual Heatmaps (Presence = 1)

if (nrow(detailed_residuals) > 0) {
  # PHASE / FASEAK ----------
  res_sig_phase <- detailed_residuals %>%
    filter(
      Test == "Phase", 
      !is.na(Category), 
      Its_presence
    ) %>%
    mutate(
      Category = trimws(toupper(Category))) %>%
    filter(Category %in% roman_levels) %>%
    mutate(
      Category = factor(
        Category, 
        levels = roman_levels, 
        ordered = TRUE)
    )
  
  if (nrow(res_sig_phase) > 0) {
    mat_res_phase <- dcast(
      res_sig_phase, 
      Category ~ Evidence, 
      value.var = "Residual", 
      fun.aggregate = mean
    )

    # Correct order by Roman numerals
    
    mat_res_phase <- mat_res_phase %>%
      mutate(Category = as.character(Category)) %>%
      arrange(
        factor(
          Category, 
          levels = roman_levels
        )
      )
    
    # Fix row names and remove the Category column 
    
    rownames(mat_res_phase) <- mat_res_phase$Category
    mat_res_phase <- mat_res_phase[
      , 
      -1, 
      drop = FALSE
    ]
    
    if (nrow(mat_res_phase) > 0) {
      png(
        file.path(
          plot_dir, 
          "Heatmap_Residuals_Presence_Phase.png"
        ), 
        width = 2400, 
        height = 1800, 
        res = 300
      )
      pheatmap(
        as.matrix(mat_res_phase),
        cluster_rows = FALSE, 
        cluster_cols = TRUE,
        show_rownames = TRUE, 
        show_colnames = TRUE,
        fontsize_row = 10, 
        fontsize_col = 9, 
        angle_col = 45,
        color = colorRampPalette(
          c(
            "blue", 
            "white", 
            "red")
        )
        (50),
        main = "Residuals (PRESENCE) Chi² per phase"
      )
      dev.off()
    }
  }
  
  
  # Region ----------------------------------------------------------------------------------------------------------------------------------------------
  
  res_sig_region <- detailed_residuals %>% filter(
    Test == "Region", 
    Its_presence
  )
  if (nrow(res_sig_region) > 0) {
    mat_res_region <- dcast(
      res_sig_region, 
      Category ~ Evidence, 
      value.var="Residual", 
      fun.aggregate=mean
    )
    rownames(mat_res_region) <- mat_res_region$Category
    mat_res_region <- mat_res_region[
      2,
      -1, 
      drop=FALSE
    ]
    png(
      file.path(
        plot_dir,
        "Heatmap_Residuals_Presence_Region.png"
      ), 
      width=2400, 
      height=1800, 
      res=300
    )
    pheatmap(
      as.matrix(mat_res_region), 
      cluster_rows = TRUE, 
      cluster_cols = TRUE,
      show_rownames = TRUE, 
      show_colnames = TRUE, 
      fontsize_row = 10, 
      fontsize_col = 9,
      angle_col = 45, 
      color = colorRampPalette(
        c(
          "blue",
          "white",
          "red"
        )
      )
      (50),
      main = "Residuals (PRESENCE) Chi² per Region"
    )
    dev.off()
  }
}
                 

# CORRESPONDENCE ANALYSIS (CA) ----------------------------------------------------------------------------------------------------------------------------------------------

tab_ca_region <- df_bin %>%
  group_by(Region) %>%
  summarise(
    across(
      all_of(evidences), 
      ~ sum(
        .==1, 
        na.rm=TRUE
      )
    ), 
    .groups="drop"
  ) %>%
  column_to_rownames("Region")
tab_ca_region <- tab_ca_region[
  rowSums(tab_ca_region)>0, 
  colSums(tab_ca_region)>0, 
  drop=FALSE
]
if (nrow(tab_ca_region)>1 && ncol(tab_ca_region)>1) {
  ca_reg <- CA(
    tab_ca_region, 
    graph=FALSE
  )
  pdf(
    file.path(
      plot_dir,
      "CA_Regions_evidences.pdf"
    ), 
    width=9, 
    height=7
  )
  print(
    fviz_ca_biplot(
      ca_reg, 
      repel=TRUE, 
      ggtheme=theme_minimal()
    )
  )
  dev.off()
}

tab_ca_phase <- df_bin %>%
  filter(!is.na(Phase_clean)) %>%
  group_by(Phase_clean) %>%
  summarise(
    across(
      all_of(evidences), 
      ~ sum(
        .==1, 
        na.rm=TRUE
      )
    ), 
    .groups="drop"
  ) %>%
  mutate(order_num = as.integer(Phase_clean)) %>%
  arrange(order_num)

rownames(tab_ca_phase) <- tab_ca_phase$Phase_clean
tab_ca_phase <- tab_ca_phase %>% dplyr::select(
  -Phase_clean, 
  -order_num
)
tab_ca_phase <- tab_ca_phase[
  rowSums(tab_ca_phase)>0, 
  colSums(tab_ca_phase)>0, 
  drop=FALSE
]
if (nrow(tab_ca_phase)>1 && ncol(tab_ca_phase)>1) {
  ca_ph <- CA(
    tab_ca_phase, 
    graph=FALSE
  )
  pdf(
    file.path(
      plot_dir,
      "CA_Phases_evidences.pdf"
    ), 
    width=9, 
    height=7
  )
  print(
    fviz_ca_biplot(
      ca_ph, 
      repel=TRUE, 
      ggtheme=theme_minimal()
    )
  )
  dev.off()
}

# MULTIPLE CORRESPONDENCE ANALYSIS (MCA) ----------------------------------------------------------------------------------------------------------------------------------------------

# Create a version of the dataset where evidences are converted to factors

mca_df <- df_bin %>%
  mutate(
    across(
      all_of(evidences), 
      ~ factor(
        ifelse(
          .==1, 
          "Yes", 
          "No")
      )
    )
  )

                 
# Select only the evidences plus grouping variables (safe and explicit)
                 
mca_df <- mca_df %>%
  dplyr::select(
    dplyr::all_of(
      c(
        evidences, 
        "Region", 
        "Phase_clean"
      )
    )
  )

const_cols <- names(
  which(
    sapply(
      mca_df[evidences], 
      nlevels
    ) < 2
  )
)
active_vars <- setdiff(
  evidences, 
  const_cols
)
if (length(active_vars) >= 2) {
  quali_sup_idx <- (length(active_vars)+1):(length(active_vars)+2)
  mca_run <- MCA(
    mca_df[
      , 
      c(
        active_vars,
        "Region",
        "Phase_clean"
      )
    ], 
    quali.sup=quali_sup_idx, 
    graph=FALSE
  )
  pdf(
    file.path(
      plot_dir,
      "MCA_evidences.pdf"
    ), 
    width=9, 
    height=7
  )
  print(
    fviz_mca_biplot(
      mca_run, 
      repel=TRUE, 
      ggtheme=theme_minimal()
    )
  )
  dev.off()
}


# PLOTS ----------------------------------------------------------------------------------------------------------------------------------------------

# Combined violin + boxplot helper
plot_violin_box <- function(
  data, 
  x_col, 
  y_col, 
  fill_pal, 
  x_lab, 
  y_lab, 
  title,
  y_scale = NULL
) {
  g <- data %>%
    filter(
      !is.na(.data[[x_col]]), 
      !is.na(.data[[y_col]])
    ) %>%
    ggplot(
      aes(
      x = .data[[x_col]], 
      y = .data[[y_col]], 
      fill = .data[[x_col]])
    ) +
    geom_violin(
      trim = FALSE, 
      color = "black", 
      alpha = 0.85
    ) +
    geom_boxplot(
      width = 0.12, 
      outlier.shape = NA, 
      alpha = 0.60
    ) +
    scale_fill_manual(
      values = fill_pal, 
      guide = "none"
    ) +
    labs(
      x = x_lab, 
      y = y_lab, 
      title = title
    ) +
    theme_pub

  if (!is.null(y_scale)) g <- g + y_scale
  g
}

difficulty_y_scale <- scale_y_continuous(
  breaks = 1:5,
  labels = c(
    "Very Low",
    "Low",ç
    "Medium",
    "Hard",
    "Very Hard"
  )
)

# Depth: violin + boxplot combined
if (!all(is.na(df_bin$Depth_m))) {
  g1 <- plot_violin_box(
    data = df_bin,
    x_col = "Region",
    y_col = "Depth_m",
    fill_pal = region_pal,
    x_lab = "Region",
    y_lab = "Reached aprox. depth (m)",
    title = "Depth by Region"
  )
  save_plot(
    "Violin_Boxplot_Depth_by_Region.png", 
    g1, 
    w = 11, 
    h = 6
  )

  g2 <- plot_violin_box(
    data = df_bin,
    x_col = "Phase_clean",
    y_col = "Depth_m",
    fill_pal = phase_pal,
    x_lab = "Phase",
    y_lab = "Reached aprox. depth (m)",
    title = "Depth by Phase"
  )
  save_plot(
    "Violin_Boxplot_Depth_by_Phase.png", 
    g2, 
    w = 9, 
    h = 6
  )
}

# Difficulty: violin + boxplot combined
if (!all(is.na(df_bin$Difficulty_num))) {
  g3 <- plot_violin_box(
    data = df_bin,
    x_col = "Region",
    y_col = "Difficulty_num",
    fill_pal = region_pal,
    x_lab = "Region",
    y_lab = "Level of difficulty of access",
    title = "Difficulty by Region",
    y_scale = difficulty_y_scale
  )
  save_plot(
    "Violin_Boxplot_Difficulty_by_Region.png", 
    g3, 
    w = 11, 
    h = 6
  )

  g4 <- plot_violin_box(
    data = df_bin,
    x_col = "Phase_clean",
    y_col = "Difficulty_num",
    fill_pal = phase_pal,
    x_lab = "Phase",
    y_lab = "Level of difficulty of access",
    title = "Difficulty by Phase",
    y_scale = difficulty_y_scale
  )
  save_plot(
    "Violin_Boxplot_Difficulty_by_Phase.png", 
    g4, 
    w = 9,
    h = 6
  )
} else if (!all(is.na(df_bin$Difficulty_fac))) {
  g5 <- df_bin %>% filter(
    !is.na(Difficulty_fac), 
    !is.na(Phase_clean)
  ) %>%
    group_by(
      Phase_clean, 
      Difficulty_fac
    ) %>% summarise(
      n=n(), 
      .groups="drop"
    ) %>%
    ggplot(
      aes(
        x=Phase_clean, 
        y=n, 
        fill=Phase_clean
      )
    ) +
    geom_col(position="fill") +
    scale_fill_manual(values = phase_pal) +
    labs(
      x="Phase", 
      y="Share", 
      fill="Phase", 
      title="Difficulty composition by Phase"
    ) +
    theme_pub
  save_plot(
    "Bars_Difficulty_by_Phase.png", 
    g5, 
    w=9, 
    h=6
  )

  g6 <- df_bin %>% filter(
    !is.na(Difficulty_fac), 
    !is.na(Region)
  ) %>%
    group_by(
      Region, 
      Difficulty_fac
    ) %>% summarise(
      n=n(), 
      .groups="drop"
    ) %>%
    ggplot(
      aes(
        x=Region, 
        y=n, 
        fill=Region
      )
    ) +
    geom_col(position="fill") +
    scale_fill_manual(values = region_pal) +
    labs(
      x="Region", 
      y="Share", 
      fill="Region", 
      title="Difficulty composition by Region"
    ) +
    theme_pub
  save_plot(
    "Bars_Difficulty_by_Region.png", 
    g6, 
    w=11, 
    h=6
  )
}


# HEATMAPS: Mean Depth by Region/Phase
                 
mat_depth_region <- df_bin %>%
  group_by(Region) %>%
  summarise(
    Mean_Depth = mean(
      Depth_m, 
      na.rm=TRUE
    ), 
    .groups="drop"
  ) %>%
  column_to_rownames("Region")
if (nrow(mat_depth_region) > 0) {
  png(
    file.path(
      plot_dir,
      "Heatmap_MeanDepth_Region.png"
    ), 
    width=2000, 
    height=1400, 
    res=300
  )
  pheatmap(
    as.matrix(mat_depth_region), 
    cluster_rows=TRUE, 
    cluster_cols=FALSE,
    main="Mean Depth by Region", 
    display_numbers=TRUE
  )
  dev.off()
}

mat_depth_phase <- df_bin %>%
  filter(!is.na(Phase_clean)) %>%
  group_by(Phase_clean) %>%
  summarise(
    Mean_Depth = mean(
      Depth_m, 
      na.rm = TRUE
    ), 
    .groups = "drop"
  ) %>%
  mutate(
    order_num = as.integer(Phase_clean)
  ) %>%
  arrange(order_num) %>%
  as.data.frame()
rownames(mat_depth_phase) <- mat_depth_phase$Phase_clean
mat_depth_phase <- mat_depth_phase %>%
  dplyr::select(
    -Phase_clean, 
    -order_num
  )
if (nrow(mat_depth_phase) > 0) {
  png(
    file.path(
      plot_dir, 
      "Heatmap_MeanDepth_Phase.png"
    ), 
    width = 2000, 
    height = 1400, 
    res = 300
  )
  pheatmap(
    mat = as.matrix(mat_depth_phase),
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    main = "Mean Depth by Phase",
    display_numbers = TRUE
  )
  dev.off()
}

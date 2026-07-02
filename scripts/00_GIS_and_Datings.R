#   SPATIAL AND CHRONOLOGICAL ANALYSIS 
# 
#   Author: Iñaki Intxaurbe Alberdi 
#   Copyright (C) 2025  Iñaki Intxaurbe
#
#   SPDX-License-Identifier: AGPL-3.0 (citation mandatory)

# Install packages and organise framework ----------------------------------------------------------------------------------------------------------------------------------------------

pkgs <- c(
  "readxl", "dplyr", "Bchron", "HDInterval","truncnorm","tidyr","purrr",
  "openxlsx","stringr","ggplot2","sf","maps","viridis","rnaturalearth",
  "rnaturalearthdata","MASS","ggridges"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(readxl)
library(dplyr)
library(Bchron)
library(HDInterval)
library(truncnorm)
library(tidyr)
library(purrr)
library(openxlsx)
library(stringr)
library(ggplot2)
library(sf)
library(maps)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(MASS)
library(ggridges)

hpd_from_draws_df <- function(
  draws, 
  cred=0.95
) 
{
  hp <- HDInterval::hdi(
    draws, 
    credMass=cred
  )
  data.frame(
    from=min(hp), 
    to=max(hp)
  )
}

parse_error_pm <- function(x) {
  if (is.na(x)) return(c(NA_real_, NA_real_))
  if (is.numeric(x)) return(c(as.numeric(x), as.numeric(x)))
  s <- gsub(
    ",",
    ".",
    as.character(x)
  )
  if (grepl("\u00B1|\\\u00B1|±", s)) {
    nums <- as.numeric(str_extract(s, "[-+]?[0-9]*\\.?[0-9]+"))
    return(c(nums, nums))
  }
  nums <- as.numeric(unlist(regmatches(s, gregexpr("[-+]?[0-9]*\\.?[0-9]+", s))))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0) return(c(NA_real_, NA_real_))
  if (length(nums) == 1) return(c(abs(nums[1]), abs(nums[1])))
  return(c(abs(nums[1]), abs(nums[2])))
}

num_safely <- function(x) {
  s <- gsub(
    ",", 
    ".", 
    as.character(x)
  )
  suppressWarnings(as.numeric(s))
}

best_phase_interval_from_pairs <- function(df_phase, max_width=10000, k=2) {
  pq_df <- df_phase %>% dplyr::filter(AQ_PQ_clean=="PQ")
  aq_df <- df_phase %>% dplyr::filter(AQ_PQ_clean=="AQ")
  if (nrow(pq_df)==0 || nrow(aq_df)==0) return(NULL)
  best <- NULL; best_w <- Inf; chosen <- NULL
  for (i in seq_len(nrow(pq_df))) {
    for (j in seq_len(nrow(aq_df))) {
      upper <- pq_df$Result[i] + k*pq_df$Error_plus[i]
      lower <- aq_df$Result[j] - k*aq_df$Error_minus[j]
      if (lower <= upper) {
        w <- upper - lower
        if (w <= max_width && w < best_w) {
          best <- c(lower, upper); best_w <- w
          chosen <- list(
            PQ_id = pq_df$`Sample ID`[i],
            AQ_id = aq_df$`Sample ID`[j],
            PQ_mean = pq_df$Result[i], PQ_plus = pq_df$Error_plus[i],
            AQ_mean = aq_df$Result[j], AQ_minus = aq_df$Error_minus[j],
            lower = lower, upper = upper, width = w
          )
        }
      }
    }
  }
  if (is.null(best)) return(NULL)
  list(interval=best, chosen=chosen)
}

round0 <- function(x) round(x, 0)

project_dir <- normalizePath(
  file.path(
    getwd(), 
    ".."
  ), 
  winslash = "/"
)

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")

dir.create(
  output_dir, 
  recursive = TRUE, 
  showWarnings = FALSE
)

file_path <- file.path(data_dir, "Table-Data-Base.xlsx")

dat <- read_excel(
  file_path, 
  sheet = "Datings"
)

dat <- dat %>%
  dplyr::filter(!is.na(Phase)) %>%
  dplyr::filter(!grepl("[?]", Phase))

dat$Result <- num_safely(dat$Result)
pm <- t(sapply(dat$Error, parse_error_pm))
dat$Error_plus  <- as.numeric(pm[,1])
dat$Error_minus <- as.numeric(pm[,2])
dat$Error_sd <- rowMeans(cbind(dat$Error_plus, dat$Error_minus), na.rm=TRUE)

dat$AQ_PQ_clean <- toupper(
    gsub(
    "\\s+", 
    "", 
    as.character(dat$`AQ/PQ`)
  )
)
dat$AQ_PQ_clean[!dat$AQ_PQ_clean %in% c("AQ","PQ","ER")] <- NA

is_c14 <- grepl(
  "C\\s*14\\s*AMS", 
  dat$Method, 
  ignore.case = TRUE
)
is_uth <- grepl(
  "U\\s*-?\\s*/?\\s*Th", 
  dat$Method, 
  ignore.case = TRUE
)


# C14 AMS ----------------------------------------------------------------------------------------------------------------------------------------------

c14 <- dat %>%
  dplyr::filter(
    is_c14, 
    !is.na(Result), 
    !is.na(Error_sd)
  ) %>%
  dplyr::mutate(
    AQ_PQ = AQ_PQ_clean,
    UniqueID = make.unique(as.character(`Sample ID`), sep = "_")
  )

bc <- BchronCalibrate(
  ages = c14$Result,
  ageSds = c14$Error_sd,
  calCurves = rep("intcal20", nrow(c14)),
  ids = c14$UniqueID,
  allowOutside = TRUE
)

apply_truncation_mask <- function(age_grid, dens, aqpq, thr) {
  if (is.na(aqpq) || aqpq == "ER") return(list(age=age_grid, dens=dens))
  if (aqpq == "AQ") mask <- age_grid >= thr
  else if (aqpq == "PQ") mask <- age_grid <= thr
  else mask <- rep(TRUE, length(age_grid))
  dens2 <- dens; dens2[!mask] <- 0
  if (sum(dens2) <= 0) dens2 <- dens
  dens2 <- dens2 / sum(dens2)
  list(age=age_grid, dens=dens2)
}

samplers_c14 <- purrr::map2(bc, seq_len(nrow(c14)), function(cal, i) {
  grid <- cal$ageGrid
  dens <- cal$densities / sum(cal$densities)
  thr  <- c14$Result[i]
  tag  <- c14$AQ_PQ[i]
  adj  <- apply_truncation_mask(grid, dens, tag, thr)
  function(n) sample(adj$age, size = n, replace = TRUE, prob = adj$dens)
})


# U/Th ----------------------------------------------------------------------------------------------------------------------------------------------

uth_raw <- dat %>%
  dplyr::filter(is_uth, !is.na(Result), !is.na(Error_plus) | !is.na(Error_minus)) %>%
  dplyr::mutate(
    AQ_PQ = AQ_PQ_clean,
    UniqueID = make.unique(as.character(`Sample ID`), sep = "_")
  )

pqaq_log <- list()
uth_bounds <- uth_raw %>% group_by(Cave, Phase) %>% group_modify(~{
  res <- best_phase_interval_from_pairs(.x, max_width=10000, k=2)
  if (is.null(res)) return(tibble())
  iv <- res$interval
  logrow <- data.frame(
    Cave = .y$Cave, Phase = .y$Phase,
    PQ_id = res$chosen$PQ_id, AQ_id = res$chosen$AQ_id,
    PQ_mean = res$chosen$PQ_mean, PQ_plus = res$chosen$PQ_plus,
    AQ_mean = res$chosen$AQ_mean, AQ_minus = res$chosen$AQ_minus,
    lower = res$chosen$lower, upper = res$chosen$upper, width = res$chosen$width
  )
  pqaq_log[[length(pqaq_log)+1]] <<- logrow
  .x %>% dplyr::mutate(Phase_min = iv[1], Phase_max = iv[2])
}) %>% ungroup()

pqaq_log_df <- if (length(pqaq_log)>0) bind_rows(pqaq_log) else data.frame()
uth <- uth_bounds

samplers_uth <- list()
if (nrow(uth) > 0) {
  samplers_uth <- lapply(seq_len(nrow(uth)), function(i) {
    mu <- uth$Result[i]
    sd <- ifelse(
      is.na(uth$Error_sd[i]), 
      mean(
        c(
          uth$Error_plus[i], 
          uth$Error_minus[i]
        ), 
      na.rm=TRUE), 
      uth$Error_sd[i]
    )
    a <- ifelse(
      !is.na(uth$Phase_min[i]), 
      uth$Phase_min[i], 
      mu - 2*uth$Error_minus[i]
    )
    b <- ifelse(
      !is.na(uth$Phase_max[i]), 
      uth$Phase_max[i], 
      mu + 2*uth$Error_plus[i]
    )
    if (!is.finite(a)) a <- mu - 4*sd
    if (!is.finite(b)) b <- mu + 4*sd
    if (a >= b) { m <- (a+b)/2; a <- m-1; b <- m+1 }
    function(n) truncnorm::rtruncnorm(n, a=a, b=b, mean=mu, sd=sd)
  })
  names(samplers_uth) <- uth$UniqueID
}


# Monte Carlo analysis -----------------------------

N <- 5000
draws_c14 <- lapply(
  samplers_c14, 
  function(f) f(N)
  ); names(draws_c14) <- c14$UniqueID
draws_uth <- lapply(
  samplers_uth, 
  function(f) f(N)
  ); names(draws_uth) <- uth$UniqueID


# Outputak / Emaitzak -------------------------------

hpd68_c14 <- lapply(draws_c14, hpd_from_draws_df, cred=0.68) %>% bind_rows()
hpd95_c14 <- lapply(draws_c14, hpd_from_draws_df, cred=0.95) %>% bind_rows()
c14_out <- c14 %>% dplyr::mutate(
  Cal_68_from = hpd68_c14$from, Cal_68_to = hpd68_c14$to,
  Cal_95_from = hpd95_c14$from, Cal_95_to = hpd95_c14$to
)

if (nrow(uth) > 0) {
  uth_out <- uth %>% dplyr::mutate(
    Cal_68_from = round0(pmax(Result - Error_minus, ifelse(is.na(Phase_min), -Inf, Phase_min))),
    Cal_68_to   = round0(pmin(Result + Error_plus,  ifelse(is.na(Phase_max),  Inf, Phase_max))),
    Cal_95_from = round0(pmax(Result - 2*Error_minus, ifelse(is.na(Phase_min), -Inf, Phase_min))),
    Cal_95_to   = round0(pmin(Result + 2*Error_plus,  ifelse(is.na(Phase_max),  Inf, Phase_max)))
  )
} else {
  uth_out <- uth
}

all_out <- bind_rows(
  c14_out %>% dplyr::mutate(Method_norm="C14 AMS"),
  uth_out %>% dplyr::mutate(Method_norm="U/Th")
)

all_draws <- c(draws_c14, draws_uth)
all_draws <- all_draws[ all_out$UniqueID ]

phases_c14 <- unique(c14_out$Phase)
phase_bounds_c14 <- lapply(phases_c14, function(ph) {
  uids <- c14_out$UniqueID[c14_out$Phase==ph]
  if (length(uids)==0) return(NULL)
  if (length(uids)==1) { d <- draws_c14[[uids]]; data.frame(phase=ph, start=d, end=d) }
  else {
    mat <- do.call(cbind, draws_c14[uids])
    start <- apply(mat,1,max); end <- apply(mat,1,min)
    data.frame(phase=ph, start=start, end=end)
  }
}) %>% bind_rows()

phase_summary_c14 <- phase_bounds_c14 %>% group_by(phase) %>% summarise(
  start_95_from = hdi(start,0.95)[2], start_95_to=hdi(start,0.95)[1],
  end_95_from = hdi(end,0.95)[2], end_95_to=hdi(end,0.95)[1],
  .groups="drop"
)

c14_out$Prob_in_Phase <- sapply(seq_len(nrow(c14_out)), function(i){
  ph <- c14_out$Phase[i]; uid <- c14_out$UniqueID[i]; d <- draws_c14[[uid]]
  pb <- phase_bounds_c14 %>% filter(phase==ph)
  mean(d >= pb$start & d <= pb$end)
})

phases_all <- unique(all_out$Phase)
phase_bounds_all <- lapply(phases_all, function(ph){
  uids <- all_out$UniqueID[all_out$Phase==ph]
  if (length(uids)==0) return(NULL)
  if (length(uids)==1) { d <- all_draws[[uids]]; data.frame(phase=ph, start=d, end=d) }
  else {
    mat <- do.call(cbind, all_draws[uids])
    start <- apply(mat,1,max); end <- apply(mat,1,min)
    data.frame(phase=ph, start=start, end=end)
  }
}) %>% bind_rows()

phase_summary_all <- phase_bounds_all %>% group_by(phase) %>% summarise(
  start_95_from = hdi(start,0.95)[2], start_95_to=hdi(start,0.95)[1],
  end_95_from = hdi(end,0.95)[2], end_95_to=hdi(end,0.95)[1],
  .groups="drop"
)

all_out$Prob_in_Phase <- sapply(seq_len(nrow(all_out)), function(i){
  ph <- all_out$Phase[i]; uid <- all_out$UniqueID[i]; d <- all_draws[[uid]]
  pb <- phase_bounds_all %>% filter(phase==ph)
  mean(d >= pb$start & d <= pb$end)
})

# Reorder phases in tables by Roman numeral / Faseak zenbakizko taula erromatarretan berrantolatzea ----
roman2int <- function(x) {
  x_clean <- gsub("Phase\\s*", "", trimws(x))
  romans <- c(
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII"
    )
  match(x_clean, romans)
}

c14_out <- c14_out %>%
  dplyr::mutate(order_num = roman2int(Phase)) %>%
  dplyr::arrange(order_num) %>%
  dplyr::select(-order_num)

if (nrow(uth) > 0) {
  uth_out <- uth_out %>%
    dplyr::mutate(order_num = roman2int(Phase)) %>%
    dplyr::arrange(order_num) %>%
    dplyr::select(-order_num)
}

all_out <- all_out %>%
  dplyr::mutate(order_num = roman2int(Phase)) %>%
  dplyr::arrange(order_num) %>%
  dplyr::select(-order_num)

phase_summary_c14 <- phase_summary_c14 %>%
  dplyr::mutate(order_num = roman2int(phase)) %>%
  dplyr::arrange(order_num) %>%
  dplyr::select(-order_num)

phase_summary_all <- phase_summary_all %>%
  dplyr::mutate(order_num = roman2int(phase)) %>%
  dplyr::arrange(order_num) %>%
  dplyr::select(-order_num)

phase_ordered <- phase_summary_all$phase


# Crearte Excel File / Excel dokumentua sortu -----------------

wb <- createWorkbook()
addWorksheet(wb, "Datings_C14_only")
addWorksheet(wb, "Phase_bounds_C14_only")
addWorksheet(wb, "Datings_combined")
addWorksheet(wb, "Phase_bounds_combined")
addWorksheet(wb, "PQ_AQ_log")

writeData(wb, "Datings_C14_only", c14_out)
writeData(wb, "Phase_bounds_C14_only", phase_summary_c14)
writeData(wb, "Datings_combined", all_out)
writeData(wb, "Phase_bounds_combined", phase_summary_all)
writeData(wb, "PQ_AQ_log", pqaq_log_df)

saveWorkbook(
  wb,
  file.path(output_dir, "Datings_results.xlsx"),
  overwrite = TRUE
  )


# Create graphs / Grafikoak sortu ------------------------

outdir <- file.path(output_dir, "Plots_results")
if (!dir.exists(outdir)) dir.create(outdir)

  # 1. Densities per phase (in Roman numerals) / Dentsitateak faseka (zenbaki erromatarretan)
  
  pdf(file.path(outdir, "Phase_densities.pdf"), width=8, height=6)
  
  for (ph in phase_ordered) {
    sub <- phase_bounds_all %>% filter(phase == ph)
    if (nrow(sub) == 0) next
    
    g1 <- ggplot(sub, aes(x = start)) +
      geom_density(fill = "skyblue", alpha = 0.5) +
      geom_vline(xintercept = HDInterval::hdi(sub$start, credMass = 0.95),
                 color = "blue", linetype = "dashed") +
      labs(title = paste("Start of phase -", ph),
           x = "Cal BP", y = "Density")
    print(g1)
    
    g2 <- ggplot(sub, aes(x = end)) +
      geom_density(fill = "salmon", alpha = 0.5) +
      geom_vline(xintercept = HDInterval::hdi(sub$end, credMass = 0.95),
                 color = "red", linetype = "dashed") +
      labs(title = paste("End of phase -", ph),
           x = "Cal BP", y = "Density")
    print(g2)
  }
  dev.off()

  # 2. Summary (in Roman numerals) / Laburpena (Zenbaki erromatarretan) 
  phase_summary_all$phase <- factor(phase_summary_all$phase, levels = phase_ordered)
  
  phase_draws_long <- phase_bounds_all %>%
    dplyr::select(phase, start, end) %>%
    tidyr::pivot_longer(cols = c(start, end),
                        names_to = "which",
                        values_to = "calBP") %>%
    dplyr::mutate(
      phase = factor(phase, levels = phase_ordered),
      which = factor(which, levels = c("start","end"))
    ) %>%
    dplyr::group_by(phase, which) %>%
    dplyr::filter(dplyr::n() > 5) %>%   
    dplyr::ungroup()
  
  plot_phase_summary <- function(phase_draws_long,
                                 phase_summary,
                                 outfile,
                                 alpha_val,
                                 scale_val,
                                 title_suffix = "") {
    
    if (nrow(phase_draws_long) < 10 ||
        length(unique(phase_draws_long$calBP)) < 2) {
      message("Skipping ridgeplot for", title_suffix,
              " (single-site or degenerate phase)")
      return(invisible(NULL))
    }
    
    p <- ggplot() +
      ggridges::geom_density_ridges(
        data = dplyr::filter(phase_draws_long, which == "start"),
        aes(x = calBP, y = phase),
        fill = "orange",
        color = "orange",
        alpha = alpha_val,
        scale = scale_val,
        rel_min_height = 0.01
      ) +
      ggridges::geom_density_ridges(
        data = dplyr::filter(phase_draws_long, which == "end"),
        aes(x = calBP, y = phase),
        fill = "#400500",
        color = "#400500",
        alpha = alpha_val,
        scale = scale_val,
        rel_min_height = 0.01
      ) +
      geom_errorbar(
        data = phase_summary,
        aes(y = phase, xmin = start_95_from, xmax = start_95_to),
        height = 0.10,
        linewidth = 0.5,
        color = "orange",
        orientation = "y"
      ) +
      geom_errorbar(
        data = phase_summary,
        aes(y = phase, xmin = end_95_from, xmax = end_95_to),
        height = 0.10,
        linewidth = 0.5,
        color = "#400500",
        orientation = "y"
      ) +
      scale_x_reverse() +
      labs(
        title = paste0("Phase posterior densities and 95% HPD intervals", title_suffix),
        x = "Cal BP",
        y = "Phase"
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
    
    ggsave(
      filename = outfile,
      plot = p,
      width = 10,
      height = 6,
      device = cairo_pdf
    )
  }
  
  # Subsets
  draws_I <- phase_draws_long %>% dplyr::filter(phase == "I")
  summary_I <- phase_summary_all %>% dplyr::filter(phase == "I")
  
  draws_rest <- phase_draws_long %>% dplyr::filter(phase != "I")
  summary_rest <- phase_summary_all %>% dplyr::filter(phase != "I")
  
  # ALL phases
  plot_phase_summary(
    phase_draws_long,
    phase_summary_all,
    alpha_val = 0.25,
    scale_val = 0.25,
    file.path(outdir, "Phase_summary_ALL.pdf")
  )
  
  # Phase I only
  plot_phase_summary(
    draws_I,
    summary_I,
    file.path(outdir, "Phase_summary_Phase_I.pdf"),
    alpha_val = 0.25,
    scale_val = 1.25,
    title_suffix = " – Phase I"
  )
  
  
  # Phases II–XII
  plot_phase_summary(
    draws_rest,
    summary_rest,
    file.path(outdir, "Phase_summary_Phase_II_to_XII.pdf"),
    alpha_val = 0.25,
    scale_val = 1.25,
    title_suffix = " – Phases II–XII"
  )
  
message("Graphs saved in the 'Plots_results' folder")


#  Spatial analyses / Analisi espazialak ---------------------
# GIS & Kernel Analysis 
  # (per-phase maps with basemap) / (fase-bakotxeko mapak basemapekin)

# Read all caves from the database / koba guztiak leiru databasetik
df_caves <- read_excel(file_path, sheet = "Caves")
df_caves$latitude  <- as.numeric(df_caves$latitude)
df_caves$longitude <- as.numeric(df_caves$longitude)
df_caves <- df_caves[!is.na(df_caves$latitude) & !is.na(df_caves$longitude), ]

# Output dir for maps / Direktorioa mapentzako

gisdir    <- file.path(
  output_dir, 
  "GIS_results"
  )
if (!dir.exists(gisdir)) dir.create(gisdir)

world_dir <- file.path(
  gisdir, 
  "Worldmaps"
  )
if (!dir.exists(world_dir)) dir.create(world_dir)

eu_dir    <- file.path(
  gisdir, 
  "European_Kernel_Densities"
  )
if (!dir.exists(eu_dir)) dir.create(eu_dir)

    # 1a) WORLD MAP (ALL CAVES, NO KERNEL) / MUNDU MAPA (KOBA DANAK; EZ KERNEL)
    
    make_world_all_caves <- function(df) {
      df$longitude <- ((df$longitude + 180) %% 360) - 180
      
      caves_ll <- st_as_sf(df, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
      world_ll <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
      world_ll <- sf::st_wrap_dateline(world_ll, options = c("WRAPDATELINE=YES"))
      
      p_world <- ggplot() +
        geom_sf(data = world_ll, fill = "gray90", color = "black", linewidth = 0.2) +
        geom_sf(data = caves_ll, color = "orange", size = 2, alpha = 0.85) +
        coord_sf(crs = "+proj=robin +lon_0=0", expand = FALSE) +   
        theme_minimal(base_size = 11) +
        labs(title = "World distribution of all caves",
             x = "Longitude", y = "Latitude")
      
      ggsave(file.path(world_dir, "World_All_Caves.png"),
             p_world, width = 12, height = 7, dpi = 300)
    }
    
    make_world_all_caves(df_caves)

    # 1b) WORLD MAP (ONLY APPROVED == YES CAVES, NO KERNEL) / MUNDU MAPA (ONARTUAK BAKARRIK == BAI KOBAK; EZ KERNEL)
    
    make_world_approved_caves <- function(df) {
      df_yes <- df %>% dplyr::filter(tolower(Approved) == "yes")
      df_yes$longitude <- ((df_yes$longitude + 180) %% 360) - 180
      
      caves_ll <- st_as_sf(df_yes, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
      world_ll <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
      world_ll <- sf::st_wrap_dateline(world_ll, options = c("WRAPDATELINE=YES"))
      
      p_world <- ggplot() +
        geom_sf(data = world_ll, fill = "gray90", color = "black", linewidth = 0.2) +
        geom_sf(data = caves_ll, color = "red", size = 2, alpha = 0.85) +
        coord_sf(crs = "+proj=robin +lon_0=0", expand = FALSE) +
        theme_minimal(base_size = 11) +
        labs(title = "World distribution of Approved caves",
             x = "Longitude", y = "Latitude")
      
      ggsave(file.path(world_dir, "World_Approved_Caves.png"),
             p_world, width = 12, height = 7, dpi = 300)
    }
    
    make_world_approved_caves(df_caves)

    # 1c) WORLD MAPS BY PHASE (APPROVED == YES & CLEAN ROMAN PHASES) / MUNDU MAPA FASEKA (ONARTURIKOAK ETA GARBITURIKO ERROMATAR ZKIAK)
    
    make_world_approved_by_phase <- function(df) {
      df_yes <- df %>% dplyr::filter(tolower(Approved) == "yes")
      
      # Filter only Approved == yes and clean phases (I, II, ... XII)
      df_yes_clean <- df_yes %>%
        dplyr::filter(grepl("^(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII)$", Phase))
      
      # Roman numerals to integers converter
            
      fases <- unique(df_yes_clean$Phase)
      fases <- fases[order(roman2int(fases))]   
      
      for (ph in fases) {
        sub <- df_yes_clean %>% dplyr::filter(Phase == ph)
        sub$longitude <- ((sub$longitude + 180) %% 360) - 180
        
        caves_ll <- st_as_sf(sub, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
        world_ll <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
        world_ll <- sf::st_wrap_dateline(world_ll, options = c("WRAPDATELINE=YES"))
        
        p_world <- ggplot() +
          geom_sf(data = world_ll, fill = "gray90", color = "black", linewidth = 0.2) +
          geom_sf(data = caves_ll, color = "red", size = 2, alpha = 0.85) +
          coord_sf(crs = "+proj=robin +lon_0=0", expand = FALSE) +
          theme_minimal(base_size = 11) +
          labs(title = paste("Approved caves - Phase", ph),
               x = "Longitude", y = "Latitude")
        
        # Numerical phase for order
        ph_num <- roman2int(ph)
        
        ggsave(file.path(world_dir, sprintf("%02d_World_Approved_Phase_%s.png", ph_num, ph)),
               p_world, width = 12, height = 7, dpi = 300)
      }
    }
    
    make_world_approved_by_phase(df_caves)

    # 2) KDE IN EUROPE (BAYESIAN CAVES ONLY, PER PHASE) / KERNEL DENTSITATE ESTIMAZIOA EUROPAN (BAYESIAN / MONTE-CARLO KOBAK BAKARRIK FASEKA)
    # Build the list of Bayesian caves (Cave+Phase present in Datings_results.xlsx)
    df_bayes <- read_excel(
     file.path(output_dir, "Datings_results.xlsx"), 
      sheet = "Datings_combined"
    ) %>%
      dplyr::distinct(Cave, Phase, .keep_all = TRUE)
    
    df_bayes_caves <- df_caves %>%
      dplyr::inner_join(df_bayes %>% dplyr::select(Cave, Phase), by = c("Cave","Phase"))
    
    make_bayes_kernels_europe <- function(df_bayes,
                                          grid_n = 280,
                                          base_gamma = 3,
                                          base_cutoff = 0.10) {
      
      # Roman numeral to integer converter / Zenbaki erromatarretik osoko bihurgailura
      roman2int <- function(x) {
        romans <- c(
          "I",
          "II",
          "III",
          "IV",
          "V",
          "VI",
          "VII",
          "VIII",
          "IX",
          "X",
          "XI",
          "XII"
          )
        match(x, romans)  
      }
      
      # Data and base map / Datuak eta base mapa
      caves_ll <- sf::st_as_sf(df_bayes, coords = c("longitude","latitude"),
                               crs = 4326, remove = FALSE)
      
      world_ll <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
      world_ll <- suppressWarnings(sf::st_make_valid(world_ll))
      world_eu <- sf::st_transform(world_ll, 3035)
      caves_eu <- sf::st_transform(caves_ll, 3035)
      
      # Robust KDE helper / KDE sendoen laguntzailea
      kde_df <- function(x, y, n = grid_n, gamma = base_gamma, cutoff = base_cutoff) {
        if (length(x) < 3L) return(NULL)
        
        rx <- range(x); ry <- range(y)
        ex <- diff(rx); ey <- diff(ry)
        
        hx <- max(ex * 0.18, 15000)
        hy <- max(ey * 0.18, 15000)
        
        if (ex < 100) x <- x + rnorm(length(x), sd = 50)
        if (ey < 100) y <- y + rnorm(length(y), sd = 50)
        
        lims <- c(rx[1] - ex*0.12, rx[2] + ex*0.12,
                  ry[1] - ey*0.12, ry[2] + ey*0.12)
        
        zz <- MASS::kde2d(x, y, h = c(hx, hy), n = n, lims = lims)
        dens <- as.vector(zz$z)
        
        rng <- range(dens, na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(NULL)
        
        dens <- (dens - rng[1]) / diff(rng)
        dens <- dens^gamma
        dens[dens < cutoff] <- NA
        
        df <- expand.grid(X = zz$x, Y = zz$y)
        df$dens <- dens
        df
      }
      
      eu_xlim <- c(2500000, 6000000)
      eu_ylim <- c(1400000, 4300000)
      
      # Loop per phase / Loop faseka
      fases <- sort(unique(df_bayes$Phase), 
                    na.last = TRUE,
                    method = "shell")
      
      for (ph in fases) {
        sub_eu <- caves_eu %>% dplyr::filter(Phase == ph)
        xy_eu  <- as.data.frame(sf::st_coordinates(sub_eu))
        if (nrow(xy_eu) > 0) names(xy_eu) <- c("X","Y")
        
        n_pts <- nrow(xy_eu)
        gamma  <- if (n_pts <= 5) 2 else base_gamma
        cutoff <- if (n_pts <= 5) 0.05 else base_cutoff
        
        kd_eu <- if (n_pts >= 3) kde_df(xy_eu$X, xy_eu$Y, gamma = gamma, cutoff = cutoff) else NULL
        
        p_eu <- ggplot() +
          geom_sf(data = world_eu, fill = "gray90", color = "black", linewidth = 0.2) +
          { if (!is.null(kd_eu))
            geom_raster(data = kd_eu, aes(X, Y, fill = dens),
                        alpha = 0.95, interpolate = TRUE) } +
          scale_fill_viridis_c(option = "magma", na.value = NA,
                               name = "Relative density") +
          geom_point(data = xy_eu, aes(X, Y), color = "red", size = 2, alpha = 0.9, na.rm = TRUE) +
          coord_sf(xlim = eu_xlim, ylim = eu_ylim, expand = FALSE, clip = "on") +
          theme_minimal(base_size = 11) +
          labs(title = paste("Bayesian Phase", ph),
               x = "X (m)", y = "Y (m)")
        
        # Roman numeral to integer for ordering / Zenbaki erromatarretatik osoetara bihurketa ordenatzeko
        ph_num <- roman2int(ph)
        if (is.na(ph_num)) ph_num <- 99 
        
        # Save files ordered according to their phase / Datuak faseen arabera gorde
        ggsave(file.path(eu_dir, sprintf("%02d_Phase_%s_Europe.png", ph_num, ph)),
               p_eu, width = 8, height = 6, dpi = 300, device = "png")
        Sys.sleep(0.1)
      }
    }
    
    make_bayes_kernels_europe(df_bayes_caves)

    # 3) KDE IN EUROPE (APPROVED == YES & CLEAN ROMAN PHASES, NOT ONLY MONTE-CARLO/BAYESIAN) 
    # 3) KDE EUROPAN (ONARTUAK ETA GARBITURIKO ERROMTARAK FASEAK, EZ BAKARRIK MONTE-CARLO/BAYESIARRAK)
    
    make_approved_kernels_europe <- function(df_caves,
                                             grid_n = 280,
                                             base_gamma = 3,
                                             base_cutoff = 0.10) {
      
      # Filter only Approved == yes and clean phases (I, II, ... XII) / Filtratu onartuak eta fase garbituak
      df_yes_clean <- df_caves %>%
        dplyr::filter(tolower(Approved) == "yes") %>%
        dplyr::filter(grepl("^(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII)$", Phase))
      
      if (nrow(df_yes_clean) == 0) {
        message("No Approved caves with clean Roman phases found")
        return(NULL)
      }
      
      # Roman numerals to integers converter / Erromatar zki-etatik osoetara bihurketa
      roman2int <- function(x) {
        romans <- c(
          "I",
          "II",
          "III",
          "IV",
          "V",
          "VI",
          "VII",
          "VIII",
          "IX",
          "X",
          "XI",
          "XII"
          )
        match(x, romans)
      }
      
      # Data and base map / Datuak eta base mapa
      caves_ll <- sf::st_as_sf(df_yes_clean, coords = c("longitude","latitude"),
                               crs = 4326, remove = FALSE)
      
      world_ll <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
      world_ll <- suppressWarnings(sf::st_make_valid(world_ll))
      world_eu <- sf::st_transform(world_ll, 3035)
      caves_eu <- sf::st_transform(caves_ll, 3035)
      
      # Helper KDE / KDE laguntzailea
      kde_df <- function(x, y, n = grid_n, gamma = base_gamma, cutoff = base_cutoff) {
        if (length(x) < 3L) return(NULL)
        rx <- range(x); ry <- range(y)
        ex <- diff(rx); ey <- diff(ry)
        hx <- max(ex * 0.18, 15000)
        hy <- max(ey * 0.18, 15000)
        if (ex < 100) x <- x + rnorm(length(x), sd = 50)
        if (ey < 100) y <- y + rnorm(length(y), sd = 50)
        lims <- c(rx[1] - ex*0.12, rx[2] + ex*0.12,
                  ry[1] - ey*0.12, ry[2] + ey*0.12)
        zz <- MASS::kde2d(x, y, h = c(hx, hy), n = n, lims = lims)
        dens <- as.vector(zz$z)
        rng <- range(dens, na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(NULL)
        dens <- (dens - rng[1]) / diff(rng)
        dens <- dens^gamma
        dens[dens < cutoff] <- NA
        df <- expand.grid(X = zz$x, Y = zz$y)
        df$dens <- dens
        df
      }
      
      eu_xlim <- c(2500000, 6000000)
      eu_ylim <- c(1400000, 4300000)
      
      # Loop per phase / Loop-ak faseka
      fases <- unique(df_yes_clean$Phase)
      fases <- fases[order(roman2int(fases))]
      
      for (ph in fases) {
        sub_eu <- caves_eu %>% dplyr::filter(Phase == ph)
        xy_eu  <- as.data.frame(sf::st_coordinates(sub_eu))
        if (nrow(xy_eu) > 0) names(xy_eu) <- c("X","Y")
        
        n_pts <- nrow(xy_eu)
        gamma  <- if (n_pts <= 5) 2 else base_gamma
        cutoff <- if (n_pts <= 5) 0.05 else base_cutoff
        
        kd_eu <- if (n_pts >= 3) kde_df(xy_eu$X, xy_eu$Y,
                                        gamma = gamma, cutoff = cutoff) else NULL
        
        p_eu <- ggplot() +
          geom_sf(data = world_eu, fill = "gray90", color = "black", linewidth = 0.2) +
          { if (!is.null(kd_eu))
            geom_raster(data = kd_eu, aes(X, Y, fill = dens),
                        alpha = 0.95, interpolate = TRUE) } +
          scale_fill_viridis_c(option = "magma", na.value = NA, name = "Relative density") +
          geom_point(data = xy_eu, aes(X, Y), color = "red", size = 2, alpha = 0.9, na.rm = TRUE) +
          coord_sf(xlim = eu_xlim, ylim = eu_ylim, expand = FALSE, clip = "on") +
          theme_minimal(base_size = 11) +
          labs(title = paste("Approved phase", ph),
               x = "X (m)", y = "Y (m)")
        
        ph_num <- roman2int(ph)
        if (is.na(ph_num)) ph_num <- 99
        
        ggsave(file.path(eu_dir, sprintf("%02d_Approved_Phase_%s_Europe.png", ph_num, ph)),
               p_eu, width = 8, height = 6, dpi = 300, device = "png")
        Sys.sleep(0.1)
      }
    }

                    
# Execute / Burutu ----------------------------------------------------

make_approved_kernels_europe(df_caves)
message("Done: 'Dating_results.xlsx' + plots in folders 'GIS_Results/' and 'Plot_Results/'.")

# OPTIONAL / AUKERAZKOA ----------------------------------------------
## Honi esker ikusi ahal diraz fase bakotxean badagozan datazino arraruak:

ph <- "X" # hemen sartu fasea

sub <- phase_bounds_all %>% filter(phase == ph)

uids <- all_out$UniqueID[all_out$Phase == ph]
mat  <- do.call(cbind, all_draws[uids])

# END ARRAROAK
idx_end <- which(sub$end < quantile(sub$end, 0.01))  # cola joven

length(idx_end)

culprit_end <- apply(mat[idx_end, , drop = FALSE], 1, function(x) {
  uids[which.min(x)]   # el más reciente manda el END
})

# START ARRAROAK

idx_start <- which(sub$start > quantile(sub$start, 0.99))  # cola vieja

length(idx_start)

culprit_start <- apply(mat[idx_start, , drop = FALSE], 1, function(x) {
  uids[which.max(x)]   # el más antiguo manda el START
})

# BALORE ARRAROAK KONTSOLAN BEGIRATU
message("BALORE ARRAROAK KONTSOLAN BEGIRATU")

table(culprit_start)
table(culprit_end)



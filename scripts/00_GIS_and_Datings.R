#   SPATIAL AND CHRONOLOGICAL ANALYSIS (CAVE-PHASE U/Th BRACKETS)
# 
#   Author: Iñaki Intxaurbe Alberdi 
#   Copyright (C) 2025  Iñaki Intxaurbe
#
#   SPDX-License-Identifier: AGPL-3.0 (citation mandatory)

# Install packages and organise framework ----------------------------------------------------------------------------------------------------------------------------------------------

pkgs <- c(
  "readxl", "dplyr", "Bchron", "HDInterval","truncnorm","tidyr","purrr",
  "openxlsx","stringr","ggplot2","sf","maps","viridis","rnaturalearth",
  "rnaturalearthdata","MASS","ggridges","coda","ArchaeoPhases","aion"
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
library(coda)
library(ArchaeoPhases)

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

base_dir <- dirname(getwd())

output_dir <- file.path(
  base_dir,
  "outputs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

tmpfile <- file.path(
  tempdir(),
  "Table_DATA.xlsx"
)

download.file(
  "https://raw.githubusercontent.com/inakiintxaurbe/Pleistocene_Human_Use_Of_Caves/refs/heads/master/data/Table-Data-Base.xlsx",
  destfile = tmpfile,
  mode = "wb"
)

dat <- readxl::read_excel(
  tmpfile,
  sheet = "Datings",
  range = readxl::cell_cols("A:K")
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

roman2int <- function(x) {
  match(x, romans)
}

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

samplers_c14 <- purrr::map(bc, function(cal) {
  grid <- cal$ageGrid
  dens <- cal$densities / sum(cal$densities)
  function(n) sample(grid, size = n, replace = TRUE, prob = dens)
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


# Monte Carlo analysis (1) -----------------------------

set.seed(666)
N <- 5000
draws_c14 <- lapply(
  samplers_c14, 
  function(f) f(N)
  ); names(draws_c14) <- c14$UniqueID
draws_uth <- lapply(
  samplers_uth, 
  function(f) f(N)
  ); names(draws_uth) <- uth$UniqueID

hpd68_c14 <- lapply(draws_c14, hpd_from_draws_df, cred=0.68) %>% bind_rows()
hpd95_c14 <- lapply(draws_c14, hpd_from_draws_df, cred=0.95) %>% bind_rows()
c14_out <- c14 %>% dplyr::mutate(
  Cal_68_from = hpd68_c14$from, Cal_68_to = hpd68_c14$to,
  Cal_95_from = hpd95_c14$from, Cal_95_to = hpd95_c14$to
)

if (nrow(uth) > 0) {
  hpd68_uth <- lapply(draws_uth, hpd_from_draws_df, cred=0.68) %>% bind_rows()
  hpd95_uth <- lapply(draws_uth, hpd_from_draws_df, cred=0.95) %>% bind_rows()
  uth_out <- uth %>% dplyr::mutate(
    Cal_68_from = round0(hpd68_uth$from),
    Cal_68_to   = round0(hpd68_uth$to),
    Cal_95_from = round0(hpd95_uth$from),
    Cal_95_to   = round0(hpd95_uth$to)
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

prob_in_phase_loo <- function(i, tab, draws) {
  uid <- tab$UniqueID[i]
  others <- tab$UniqueID[tab$Phase == tab$Phase[i] & tab$UniqueID != uid]
  if (length(others) == 0) return(NA_real_)
  other_mat <- do.call(cbind, draws[others])
  if (is.null(dim(other_mat))) other_mat <- matrix(other_mat, ncol=1)
  lower <- apply(other_mat, 1, min)
  upper <- apply(other_mat, 1, max)
  d <- draws[[uid]]
  mean(d >= lower & d <= upper)
}

phase_summary_c14 <- phase_bounds_c14 %>% group_by(phase) %>% summarise(
  start_95_from = hdi(start,0.95)[1], start_95_to=hdi(start,0.95)[2],
  end_95_from = hdi(end,0.95)[1], end_95_to=hdi(end,0.95)[2],
  .groups="drop"
)

c14_out$Prob_in_Phase <- vapply(
  seq_len(nrow(c14_out)), prob_in_phase_loo, numeric(1),
  tab = c14_out, draws = draws_c14
)


# Monte Carlo analysis (2) -----------------------------

uth_bracket_draws <- list()
uth_bracket_log <- list()
if (nrow(pqaq_log_df) > 0) {
  for (j in seq_len(nrow(pqaq_log_df))) {
    z <- pqaq_log_df[j, , drop=FALSE]
    key <- paste(z$Cave, z$Phase, sep="\r")
    pq_sd <- uth_raw$Error_sd[which(uth_raw$Cave==z$Cave &
      uth_raw$Phase==z$Phase & uth_raw$`Sample ID`==z$PQ_id &
      uth_raw$AQ_PQ=="PQ")[1]]
    aq_sd <- uth_raw$Error_sd[which(uth_raw$Cave==z$Cave &
      uth_raw$Phase==z$Phase & uth_raw$`Sample ID`==z$AQ_id &
      uth_raw$AQ_PQ=="AQ")[1]]
    if (!is.finite(pq_sd) || pq_sd <= 0) pq_sd <- 1
    if (!is.finite(aq_sd) || aq_sd <= 0) aq_sd <- 1
    upper <- stats::rnorm(N, mean=z$PQ_mean, sd=pq_sd)
    lower <- stats::rnorm(N, mean=z$AQ_mean, sd=aq_sd)
    valid <- is.finite(lower) & is.finite(upper) & lower <= upper
    uth_bracket_draws[[key]] <- list(lower=lower, upper=upper, valid=valid)
    uth_bracket_log[[length(uth_bracket_log)+1]] <- data.frame(
      Cave=z$Cave, Phase=z$Phase, N_PQ=1L, N_AQ=1L,
      Bracket_consistency=mean(valid),
      AQ_younger_median=stats::median(lower),
      PQ_older_median=stats::median(upper),
      stringsAsFactors=FALSE
    )
  }
}
uth_bracket_log_df <- if (length(uth_bracket_log)) dplyr::bind_rows(uth_bracket_log) else
  data.frame(
    Cave=character(),
    Phase=character(),
    N_PQ=integer(),
    N_AQ=integer(),
    Bracket_consistency=numeric(),
    AQ_younger_median=numeric(),
    PQ_older_median=numeric()
    )


# Chronologies / Kronologiak ----------------------------------------

cave_phase_keys <- unique(all_out[,c("Cave","Phase")])
cave_phase_rows <- list()
phase_event_draws <- draws_c14
bounded_event_draws <- list()
for (j in seq_len(nrow(cave_phase_keys))) {
  cave <- as.character(cave_phase_keys$Cave[j]); ph <- as.character(cave_phase_keys$Phase[j])
  key <- paste(cave,ph,sep="\r")
  ids <- c14$UniqueID[c14$Cave==cave & c14$Phase==ph]
  bracket <- uth_bracket_draws[[key]]
  n_c14 <- length(ids)
  has_bracket <- !is.null(bracket)
  bracket_ok <- if (has_bracket) mean(bracket$valid) else NA_real_
  compatibility <- NA_real_
  if (n_c14 && has_bracket) {
    compatibility <- mean(vapply(ids, function(uid) {
      x <- draws_c14[[uid]]
      mean(bracket$valid & x >= bracket$lower & x <= bracket$upper)
    }, numeric(1)))
  }
  event <- NULL
  event_kind <- "No direct event date"
  if (n_c14) {
    event_kind <- "Direct C14"
    m <- do.call(cbind, draws_c14[ids])
    if (is.null(dim(m))) m <- matrix(m,ncol=1)
    event <- rowMeans(m) # descriptive central event, not a phase boundary
  } else if (has_bracket && any(bracket$valid)) {
    event_kind <- "UTh bracket proxy (not direct)"
    ok <- which(bracket$valid)
    chosen <- sample(ok, size=N, replace=TRUE)
    event <- stats::runif(N, min=bracket$lower[chosen], max=bracket$upper[chosen])
    uid <- paste0("UTH_BOUND_EVENT_", j)
    bounded_event_draws[[uid]] <- event
    phase_event_draws[[uid]] <- event
  }
  if (!is.null(event)) {
    ev68 <- hpd_from_draws_df(event,cred=.68)
    ev95 <- hpd_from_draws_df(event,cred=.95)
    event_median <- stats::median(event)
  } else {
    ev68 <- data.frame(from=NA_real_,to=NA_real_)
    ev95 <- ev68
    event_median <- NA_real_
  }
  cave_phase_rows[[j]] <- data.frame(
    Cave=cave, Phase=ph, N_C14=n_c14,
    N_UTh=sum(uth_raw$Cave==cave & uth_raw$Phase==ph),
    Has_AQ_PQ_bracket=has_bracket,
    Bracket_consistency=bracket_ok,
    C14_within_UTh_bracket=compatibility,
    Chronological_basis=event_kind,
    Event_median_calBP=event_median,
    Event_68_from=ev68$from, Event_68_to=ev68$to,
    Event_95_from=ev95$from, Event_95_to=ev95$to,
    stringsAsFactors=FALSE
  )
}
cave_phase_chronology <- dplyr::bind_rows(cave_phase_rows) %>%
  dplyr::mutate(order_num=match(Phase,romans)) %>%
  dplyr::arrange(order_num,Cave) %>% dplyr::select(-order_num)


# Envelopes / Gainjarketak ---------------------------------------

phase_bounds_all <- lapply(unique(cave_phase_chronology$Phase), function(ph) {
  ids <- c14$UniqueID[c14$Phase==ph]
  proxy_keys <- which(cave_phase_chronology$Phase==ph &
                      cave_phase_chronology$Chronological_basis=="UTh bracket proxy (not direct)")
  proxy_ids <- character()
  if (length(proxy_keys)) {
    for (j in proxy_keys) {
      row <- which(cave_phase_keys$Cave==cave_phase_chronology$Cave[j] &
                   cave_phase_keys$Phase==ph)
      proxy_ids <- c(proxy_ids,paste0("UTH_BOUND_EVENT_",row[1]))
    }
  }
  all_ids <- c(ids,proxy_ids)
  all_ids <- all_ids[all_ids %in% names(phase_event_draws)]
  if (!length(all_ids)) return(NULL)
  m <- do.call(cbind, phase_event_draws[all_ids])
  if (is.null(dim(m))) m <- matrix(m,ncol=1)
  data.frame(phase=ph,start=apply(m,1,max),end=apply(m,1,min))
}) %>% dplyr::bind_rows()

phase_summary_all <- phase_bounds_all %>% dplyr::group_by(phase) %>% dplyr::summarise(
  start_95_from=HDInterval::hdi(start,.95)[1], start_95_to=HDInterval::hdi(start,.95)[2],
  end_95_from=HDInterval::hdi(end,.95)[1], end_95_to=HDInterval::hdi(end,.95)[2],
  .groups="drop"
)

all_out$Prob_in_Phase <- NA_real_
c14_indices <- which(all_out$Method_norm == "C14 AMS")
if (length(c14_indices)) {
  all_out$Prob_in_Phase[c14_indices] <- vapply(c14_indices, function(i) {
    uid <- all_out$UniqueID[i]
    others <- c14$UniqueID[c14$Phase == all_out$Phase[i] & c14$UniqueID != uid]
    if (!length(others)) return(NA_real_)
    m <- do.call(cbind,draws_c14[others])
    if (is.null(dim(m))) m <- matrix(m,ncol=1)
    d <- draws_c14[[uid]]
    mean(d >= apply(m,1,min) & d <= apply(m,1,max))
  },numeric(1))
}

phase_overlap <- list()
for (j in seq_len(length(romans)-1)) {
  a <- romans[j]; b <- romans[j+1]
  aa <- phase_bounds_all[phase_bounds_all$phase==a, ,drop=FALSE]
  bb <- phase_bounds_all[phase_bounds_all$phase==b, ,drop=FALSE]
  if (!nrow(aa) || !nrow(bb)) next
  overlap <- pmax(0,pmin(aa$start,bb$start)-pmax(aa$end,bb$end))
  gap <- pmax(0,pmax(aa$end-bb$start,bb$end-aa$start))
  phase_overlap[[length(phase_overlap)+1]] <- data.frame(
    Phase_older=a,Phase_younger=b,
    P_overlap=mean(overlap>0),
    Overlap_median_calYears=median(overlap),
    Gap_median_calYears=median(gap),
    P_nonoverlap=mean(overlap==0),
    stringsAsFactors=FALSE
  )
}
phase_overlap_df <- if (length(phase_overlap)) dplyr::bind_rows(phase_overlap) else
  data.frame(Phase_older=character(),Phase_younger=character(),P_overlap=numeric(),
             Overlap_median_calYears=numeric(),Gap_median_calYears=numeric(),P_nonoverlap=numeric())

chronology_methods <- data.frame(
  Item=c(
    "Inference",
    "UTh AQ/PQ",
    "Cave/phase",
    "Phase envelopes",
    "Overlap",
    "MCMC diagnostics"
    ),
  Interpretation=c(
    "Original independent Monte Carlo sheets; separate hierarchical MCMC sheets",
    "Within-cave/phase bounds on events, not direct phase beginning/end dates",
    "Direct C14 when available; otherwise exploratory uniform UTh bracket event proxy",
    "Sampled earliest/latest dated events, not posterior phase start/end parameters",
    "Descriptive sampled-envelope overlap, not proof of occupation continuity or absence",
    "See MCMC_Diagnostics for four-chain Metropolis-Hastings diagnostics"
  ),stringsAsFactors=FALSE
)


# Bayesian MCMC ------------------------------------------------

mcmc_chains <- 4L
mcmc_burn <- 2500L
mcmc_keep <- 3000L
mcmc_thin <- 2L
mcmc_phase_sd <- 3000 
mcmc_step_event <- 350
mcmc_step_center <- 500

mcmc_events <- c14 %>% dplyr::transmute(
  EventID=UniqueID, Cave=as.character(Cave), Phase=as.character(Phase),
  Kind="C14", SourceID=UniqueID
)
if (nrow(uth_bracket_log_df)) {
  mcmc_events <- dplyr::bind_rows(mcmc_events,
    uth_bracket_log_df %>% dplyr::transmute(
      EventID=paste0("UTH_EVENT_",seq_len(dplyr::n())),
      Cave=as.character(Cave), Phase=as.character(Phase),
      Kind="UTh_bracket", SourceID=EventID
    ))
}

mcmc_events$Phase <- as.character(mcmc_events$Phase)
mcmc_phases <- unique(mcmc_events$Phase)
mcmc_phases <- romans[romans %in% mcmc_phases]
mcmc_events$phase_index <- match(mcmc_events$Phase,mcmc_phases)
mcmc_age_max <- max(c(70000, c14$Result+10*c14$Error_sd,
                            uth_raw$Result+10*uth_raw$Error_sd),na.rm=TRUE)
mcmc_age_max <- ceiling(mcmc_age_max/1000)*1000

# U/Th likelihood uses only the selected narrowest valid AQ/PQ pair.
mcmc_loglik <- vector("list",nrow(mcmc_events))
for (i in seq_len(nrow(mcmc_events))) {
  if (mcmc_events$Kind[i]=="C14") {
    cal <- bc[[match(mcmc_events$SourceID[i],c14$UniqueID)]]
    xx <- as.numeric(cal$ageGrid)
    yy <- as.numeric(cal$densities)
    good <- is.finite(xx)&is.finite(yy)&yy>0
    xx <- xx[good]; yy <- yy[good]
    oo <- order(xx); xx <- xx[oo]; yy <- yy[oo]
    mcmc_loglik[[i]] <- local({
      x <- xx; y <- yy
      function(t) {
        if (!is.finite(t)||t<min(x)||t>max(x)) return(-Inf)
        d <- stats::approx(x,y,xout=t,rule=1)$y
        if (!is.finite(d)||d<=0) return(-Inf)
        log(d)
      }
    })
  } else {
    z <- pqaq_log_df[as.character(pqaq_log_df$Cave)==mcmc_events$Cave[i] &
                     as.character(pqaq_log_df$Phase)==mcmc_events$Phase[i],,drop=FALSE]
    mcmc_loglik[[i]] <- local({
      pqmu <- z$PQ_mean[1]; aqmu <- z$AQ_mean[1]
      pqsd <- uth_raw$Error_sd[which(uth_raw$Cave==z$Cave[1] &
        uth_raw$Phase==z$Phase[1] & uth_raw$`Sample ID`==z$PQ_id[1] &
        uth_raw$AQ_PQ=="PQ")[1]]
      aqsd <- uth_raw$Error_sd[which(uth_raw$Cave==z$Cave[1] &
        uth_raw$Phase==z$Phase[1] & uth_raw$`Sample ID`==z$AQ_id[1] &
        uth_raw$AQ_PQ=="AQ")[1]]
      if (!is.finite(pqsd) || pqsd <= 0) pqsd <- 1
      if (!is.finite(aqsd) || aqsd <= 0) aqsd <- 1
      function(t) {
        if (!is.finite(t)) return(-Inf)
        stats::pnorm(t,mean=pqmu,sd=pqsd,lower.tail=TRUE,log.p=TRUE) +
          stats::pnorm(t,mean=aqmu,sd=aqsd,lower.tail=FALSE,log.p=TRUE)
      }
    })
  }
}

mcmc_initial <- function(i) {
  if (mcmc_events$Kind[i]=="C14") {
    cal <- bc[[match(mcmc_events$SourceID[i],c14$UniqueID)]]
    return(sample(cal$ageGrid,1,prob=cal$densities))
  }
  z <- pqaq_log_df[as.character(pqaq_log_df$Cave)==mcmc_events$Cave[i] &
                   as.character(pqaq_log_df$Phase)==mcmc_events$Phase[i],,drop=FALSE]
  (z$PQ_mean[1]+z$AQ_mean[1])/2
}

mcmc_event_names <- make.names(mcmc_events$EventID,unique=TRUE)
mcmc_phase_names <- paste0("Center_",mcmc_phases)
mcmc_draws <- vector("list",mcmc_chains)
mcmc_accept <- matrix(0,mcmc_chains,2,dimnames=list(NULL,c("Events","Centers")))
for (ch in seq_len(mcmc_chains)) {
  set.seed(2026L+ch)
  age <- vapply(seq_len(nrow(mcmc_events)),mcmc_initial,numeric(1))
  center <- vapply(seq_along(mcmc_phases),function(j)
    stats::median(age[mcmc_events$phase_index==j]),numeric(1))
  if (ch>1L) {
    age <- pmax(0,pmin(mcmc_age_max,age+stats::rnorm(length(age),0,300)))
    center <- pmax(0,pmin(mcmc_age_max,center+stats::rnorm(length(center),0,500)))
  }
  # Ensure all chains start on finite support of the calibrated likelihood.
  for (i in seq_along(age)) if (!is.finite(mcmc_loglik[[i]](age[i]))) age[i] <- mcmc_initial(i)
  total <- mcmc_burn + mcmc_keep*mcmc_thin
  out <- matrix(NA_real_,nrow=mcmc_keep,
                ncol=length(age)+length(center),
                dimnames=list(NULL,c(mcmc_event_names,mcmc_phase_names)))
  saved <- 0L
  for (iter in seq_len(total)) {
    for (i in seq_along(age)) {
      p <- mcmc_events$phase_index[i]
      candidate <- age[i]+stats::rnorm(1,0,mcmc_step_event)
      if (candidate<0||candidate>mcmc_age_max) next
      logratio <- mcmc_loglik[[i]](candidate)-mcmc_loglik[[i]](age[i])+
        stats::dnorm(candidate,center[p],mcmc_phase_sd,log=TRUE)-
        stats::dnorm(age[i],center[p],mcmc_phase_sd,log=TRUE)
      if (is.finite(logratio) && log(stats::runif(1))<logratio) {
        age[i] <- candidate
        mcmc_accept[ch,1] <- mcmc_accept[ch,1]+1
      }
    }
    for (j in seq_along(center)) {
      candidate <- center[j]+stats::rnorm(1,0,mcmc_step_center)
      if (candidate<0||candidate>mcmc_age_max) next
      ix <- which(mcmc_events$phase_index==j)
      logratio <- sum(stats::dnorm(age[ix],candidate,mcmc_phase_sd,log=TRUE)-
                      stats::dnorm(age[ix],center[j],mcmc_phase_sd,log=TRUE))
      if (is.finite(logratio) && log(stats::runif(1))<logratio) {
        center[j] <- candidate
        mcmc_accept[ch,2] <- mcmc_accept[ch,2]+1
      }
    }
    if (iter>mcmc_burn && (iter-mcmc_burn)%%mcmc_thin==0L) {
      saved <- saved+1L
      out[saved,] <- c(age,center)
    }
  }
  mcmc_draws[[ch]] <- out
}
mcmc_coda <- coda::mcmc.list(lapply(mcmc_draws,coda::mcmc))
mcmc_ess <- coda::effectiveSize(mcmc_coda)
mcmc_rhat <- tryCatch(coda::gelman.diag(mcmc_coda,multivariate=FALSE,
                                        autoburnin=FALSE)$psrf[,1],
                       error=function(e) rep(NA_real_,ncol(mcmc_draws[[1]])))
mcmc_diagnostics <- data.frame(
  Parameter=colnames(mcmc_draws[[1]]),
  ESS=as.numeric(mcmc_ess[colnames(mcmc_draws[[1]])]),
  Rhat=as.numeric(mcmc_rhat[colnames(mcmc_draws[[1]])]),
  stringsAsFactors=FALSE
)
mcmc_diagnostics$Converged <- is.finite(mcmc_diagnostics$Rhat) &
  mcmc_diagnostics$Rhat<1.05 & mcmc_diagnostics$ESS>=400
mcmc_pooled <- do.call(rbind,mcmc_draws)

mcmc_cave_phase <- lapply(seq_len(nrow(cave_phase_chronology)),function(i) {
  ix <- which(mcmc_events$Cave==as.character(cave_phase_chronology$Cave[i]) &
              mcmc_events$Phase==as.character(cave_phase_chronology$Phase[i]))
  if (!length(ix)) return(NULL)
  mm <- mcmc_pooled[,ix,drop=FALSE]
  oldest <- apply(mm,1,max); youngest <- apply(mm,1,min)
  data.frame(Cave=cave_phase_chronology$Cave[i],Phase=cave_phase_chronology$Phase[i],
    N_Events=length(ix),
    Oldest_95_from=unname(stats::quantile(oldest,.025)),
    Oldest_95_to=unname(stats::quantile(oldest,.975)),
    Youngest_95_from=unname(stats::quantile(youngest,.025)),
    Youngest_95_to=unname(stats::quantile(youngest,.975)),
    Diagnostics_OK=all(mcmc_diagnostics$Converged[ix]),stringsAsFactors=FALSE)
}) %>% dplyr::bind_rows()

mcmc_phase_ranges <- lapply(mcmc_phases,function(ph) {
  ix <- which(mcmc_events$Phase==ph)
  mm <- mcmc_pooled[,ix,drop=FALSE]
  data.frame(phase=ph,start=apply(mm,1,max),end=apply(mm,1,min))
})
names(mcmc_phase_ranges) <- mcmc_phases
mcmc_phase_summary <- dplyr::bind_rows(lapply(mcmc_phases,function(ph) {
  x <- mcmc_phase_ranges[[ph]]
  data.frame(phase=ph,
    start_95_from=unname(stats::quantile(x$start,.025)),
    start_95_to=unname(stats::quantile(x$start,.975)),
    end_95_from=unname(stats::quantile(x$end,.025)),
    end_95_to=unname(stats::quantile(x$end,.975)),
    Diagnostics_OK=all(mcmc_diagnostics$Converged[c(which(mcmc_events$Phase==ph),
                         length(mcmc_event_names)+match(ph,mcmc_phases))]),
    stringsAsFactors=FALSE)
}))

mcmc_overlap <- list()
for (j in seq_len(length(romans)-1L)) {
  a <- romans[j]; b <- romans[j+1L]
  if (!a %in% mcmc_phases || !b %in% mcmc_phases) next
  aa <- mcmc_phase_ranges[[a]]; bb <- mcmc_phase_ranges[[b]]
  overlap <- pmax(0,pmin(aa$start,bb$start)-pmax(aa$end,bb$end))
  gap <- pmax(0,pmax(aa$end-bb$start,bb$end-aa$start))
  mcmc_overlap[[length(mcmc_overlap)+1L]] <- data.frame(
    Phase_older=a,Phase_younger=b,P_overlap=mean(overlap>0),
    Overlap_median_calYears=stats::median(overlap),
    Gap_median_calYears=stats::median(gap),
    Gap_if_present_median_calYears=if(any(gap>0)) stats::median(gap[gap>0]) else NA_real_,
    P_nonoverlap=mean(gap>0),stringsAsFactors=FALSE)
}
mcmc_phase_overlap <- dplyr::bind_rows(mcmc_overlap)

# ArchaeoPhases 2.x: post-process the actual MCMC phase-range draws.
mcmc_phase_matrix <- do.call(cbind,lapply(mcmc_phases,function(ph)
  as.matrix(mcmc_phase_ranges[[ph]][,c("start","end")])))
colnames(mcmc_phase_matrix) <- as.vector(rbind(paste0(mcmc_phases,"_start"),
                                                paste0(mcmc_phases,"_end")))
mcmc_archaeophases <- ArchaeoPhases::as_phases(
  mcmc_phase_matrix,calendar=aion::BP(),
  start=seq(1L,ncol(mcmc_phase_matrix),by=2L),
  stop=seq(2L,ncol(mcmc_phase_matrix),by=2L),names=mcmc_phases)
mcmc_archaeophases_boundaries <- as.data.frame(
  ArchaeoPhases::boundaries(mcmc_archaeophases,level=.95))
mcmc_acceptance <- data.frame(Chain=seq_len(mcmc_chains),
  Events=mcmc_accept[,1]/(mcmc_burn+mcmc_keep*mcmc_thin)/nrow(mcmc_events),
  Centers=mcmc_accept[,2]/(mcmc_burn+mcmc_keep*mcmc_thin)/length(mcmc_phases))


# Output-ak -------------------------------------------------------------

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
addWorksheet(wb, "Cave_Phase_Chronology")
addWorksheet(wb, "UTh_Brackets")
addWorksheet(wb, "Phase_Overlap")
addWorksheet(wb, "Chronology_Methods")
addWorksheet(wb, "MCMC_Cave_Phase")
addWorksheet(wb, "MCMC_Phase_bounds")
addWorksheet(wb, "MCMC_Phase_Overlap")
addWorksheet(wb, "MCMC_Diagnostics")
addWorksheet(wb, "MCMC_Acceptance")
addWorksheet(wb, "ArchaeoPhases_Boundaries")

writeData(wb, "Datings_C14_only", c14_out)
writeData(wb, "Phase_bounds_C14_only", phase_summary_c14)
writeData(wb, "Datings_combined", all_out)
writeData(wb, "Phase_bounds_combined", phase_summary_all)
writeData(wb, "PQ_AQ_log", pqaq_log_df)
writeData(wb, "Cave_Phase_Chronology", cave_phase_chronology)
writeData(wb, "UTh_Brackets", uth_bracket_log_df)
writeData(wb, "Phase_Overlap", phase_overlap_df)
writeData(wb, "Chronology_Methods", chronology_methods)
writeData(wb, "MCMC_Cave_Phase", mcmc_cave_phase)
writeData(wb, "MCMC_Phase_bounds", mcmc_phase_summary)
writeData(wb, "MCMC_Phase_Overlap", mcmc_phase_overlap)
writeData(wb, "MCMC_Diagnostics", mcmc_diagnostics)
writeData(wb, "MCMC_Acceptance", mcmc_acceptance)
writeData(wb, "ArchaeoPhases_Boundaries", mcmc_archaeophases_boundaries)

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
df_caves <- read_excel(tmpfile, sheet = "Caves")
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

table_1 <- table(culprit_start)
table_2 <- table(culprit_end)

print(table_1)
print(table_2)

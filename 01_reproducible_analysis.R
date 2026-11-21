# =============================================================================
# Reproducible analysis code
# Risk of bias in psychotherapy trials for depression: a large-scale
# empirical examination using a structured assessment tool developed for
# psychotherapy research
#
# This single script reproduces every statistical result reported in the
# manuscript (Objectives 1-5) from the two shared data files:
#   - irr_ratings.csv     (double-rating subsample, 24 RoB items)
#   - analysis_data.csv   (941 instrument-level rows: RoB ratings,
#                                   effect sizes, and the trial-characteristic
#                                   predictors used in Objectives 2-5)
#
# =============================================================================

library(readxl)   # unused here, kept only if you swap the CSVs for .xlsx
library(dplyr)
library(tidyr)
library(MASS)       # polr() - ordinal logistic regression
library(psych)       # cohen.kappa() - CI for Objective 1 kappa
library(metafor)     # rma.mv() - CHE meta-analytic models
library(clubSandwich) # CR2 robust variance estimation

select <- dplyr::select
filter <- dplyr::filter

data_dir <- "shared_data"

rob_levels <- c("low risk", "some concerns", "high risk")

domain_labels <- c(
  rob_study_lvl = "Overall RoB score",
  rob_d1 = "Domain 1 - Randomization process",
  rob_d2 = "Domain 2 - Deviations from intended interventions",
  rob_d3 = "Domain 3 - Missing outcome data",
  rob_d4 = "Domain 4 - Measurement of the outcome",
  rob_d5 = "Domain 5 - Selective outcome reporting"
)

# =============================================================================
# OBJECTIVE 1 - Inter-rater reliability
# =============================================================================

irr <- read.csv(file.path(data_dir, "irr_ratings.csv"), na.strings = "")

rating_cols <- c(
  "d1_1","d1_2","d1_3","d1_4",
  "d2_5","d2_6","d2_7","d2_8","d2_9",
  "d3_10","d3_11","d3_12","d3_13","d3_14",
  "d4_15","d4_16","d4_17","d4_18",
  "d5_19","d5_21","d5_22","d5_23","d5_24"   # d5_20 excluded: free-text field
)

# Sensitivity analysis: the 16 manually-rated items only, excluding the 7
# items that are pre-filled automatically from already-extracted data
# (d1_3, d1_4, d3_10, d3_12, d3_13, d4_16, d4_17) and therefore identical by
# construction for both raters.
manual_cols <- c(
  "d1_1", "d1_2",
  "d2_5", "d2_6", "d2_7", "d2_8", "d2_9",
  "d3_11", "d3_14",
  "d4_15", "d4_18",
  "d5_19", "d5_21", "d5_22", "d5_23", "d5_24"
)

pairs <- list(
  list(group = "CM & TF",  rater1 = "CM", rater2 = "TF"),
  list(group = "CM & EsK", rater1 = "CM", rater2 = "EsK"),
  list(group = "AS & PK",  rater1 = "AS", rater2 = "PK"),
  list(group = "EK & OP",  rater1 = "EK", rater2 = "OP"),
  list(group = "ES & MS",  rater1 = "ES", rater2 = "MS"),
  list(group = "LT & EsK", rater1 = "LT", rater2 = "EsK"),
  list(group = "PC & YW",  rater1 = "PC", rater2 = "YW"),
  list(group = "RT & AT",  rater1 = "RT", rater2 = "AT"),
  list(group = "RT & KY",  rater1 = "RT", rater2 = "KY"),
  list(group = "CM & MC",  rater1 = "CM", rater2 = "MC")
)

get_rater <- function(grp, id) {
  irr %>% filter(trimws(group) == grp, trimws(rater) == id)
}

# Binarised item-level Cohen's kappa + SE (Yes/PY vs. No/PN or NI)
item_kappa <- function(d1, d2, cols) {
  r1 <- as.matrix(d1[, cols]); r2 <- as.matrix(d2[, cols])
  r1[!r1 %in% c("Yes/PY", "No/PN", "NI")] <- NA
  r2[!r2 %in% c("Yes/PY", "No/PN", "NI")] <- NA
  keep <- !is.na(r1) & !is.na(r2)
  r1b <- (r1[keep] == "Yes/PY"); r2b <- (r2[keep] == "Yes/PY")
  obs_agree <- round(mean(r1b == r2b) * 100, 1)
  se <- tryCatch({
    ck <- psych::cohen.kappa(cbind(factor(r1b, levels = c(FALSE, TRUE)),
                                    factor(r2b, levels = c(FALSE, TRUE))))
    (ck$confid["unweighted kappa", "upper"] - ck$confid["unweighted kappa", "lower"]) / (2 * 1.96)
  }, error = function(e) NA_real_)
  p0 <- mean(r1b == r2b); p1 <- mean(r1b); p2 <- mean(r2b)
  pe <- p1 * p2 + (1 - p1) * (1 - p2)
  list(kappa = (p0 - pe) / (1 - pe), se = se, obs_agree = obs_agree)
}

# Score-level agreement: Gwet's AC2, ordinal weights (requires irrCAC)
score_ac2 <- function(r1, r2) {
  r1 <- tolower(trimws(r1)); r2 <- tolower(trimws(r2))
  valid <- r1 %in% rob_levels & r2 %in% rob_levels
  if (sum(valid) < 2) return(list(ac2 = NA_real_, obs_agree = NA_real_))
  obs_agree <- round(mean(r1[valid] == r2[valid]) * 100, 1)
  ac2 <- tryCatch(
    irrCAC::gwet.ac1.raw(data.frame(r1 = r1[valid], r2 = r2[valid]),
                          weights = "ordinal", categ.labels = rob_levels)$est$coeff.val,
    error = function(e) NA_real_
  )
  list(ac2 = ac2, obs_agree = obs_agree)
}

# Study-level RoB scores are derived from the 24 item ratings using the
# same scoring algorithms as the rest of the manuscript (metapsyTools, cited
# in Methods/References). This calls the package rather than re-deriving the
# domain/overall scoring rules here, since those rules live in Supplement B
# and must stay in sync with the single implementation the whole project
# uses - duplicating them in this script would risk the two silently
# diverging. Install metapsyTools (see manuscript references) before running
# this section.
score_from_items <- function(d) {
  # createRobRatings() expects the five *_notes columns to be present, even
  # though their free-text content is never used in scoring (which is why
  # they were dropped from the shared data file). Add them back as empty
  # placeholders so the function runs unchanged.
  for (nc in c("d1_notes", "d2_notes", "d3_notes", "d4_notes", "d5_notes", "d5_20")) {
    if (!nc %in% names(d)) d[[nc]] <- NA_character_
  }
  metapsyTools:::createRobRatings(rob.data = d, database = NULL)$rob.data
}

# Kappa and AC2 are both bounded in [-1, 1] (1 = perfect agreement). A plain
# symmetric t-interval doesn't know that, and can overshoot the boundary when
# several pairs score at or near it (e.g. several pairs at kappa/AC2 = 1.000
# on a low-disagreement domain) - clip to the coefficient's valid range.
naive_ci <- function(x, bounds = c(-1, 1)) {
  x <- x[!is.na(x)]; k <- length(x)
  se <- sd(x) / sqrt(k)
  lower <- max(mean(x) - qt(.975, k - 1) * se, bounds[1])
  upper <- min(mean(x) + qt(.975, k - 1) * se, bounds[2])
  c(mean = mean(x), lower = lower, upper = upper)
}

cat("=== OBJECTIVE 1: Inter-rater reliability ===\n")
tryCatch({
  irr_item_rows  <- list()
  irr_score_rows <- list()

  for (p in pairs) {
    d1 <- get_rater(p$group, p$rater1); d2 <- get_rater(p$group, p$rater2)
    if (nrow(d1) == 0 || nrow(d2) == 0) next
    ik     <- item_kappa(d1, d2, rating_cols)
    ik_man <- item_kappa(d1, d2, manual_cols)
    irr_item_rows[[p$group]] <- data.frame(pair = p$group, kappa = ik$kappa,
                                            obs_agree = ik$obs_agree,
                                            kappa_manual = ik_man$kappa,
                                            obs_agree_manual = ik_man$obs_agree)

    score_cols <- c("d1", "d2", "d3", "d4", "d5", "rob")
    s1 <- score_from_items(d1)[, c("study", score_cols)]
    s2 <- score_from_items(d2)[, c("study", score_cols)]
    m  <- merge(s1, s2, by = "study", suffixes = c("_r1", "_r2"))
    row <- data.frame(pair = p$group)
    for (col in score_cols) {
      r <- score_ac2(m[[paste0(col, "_r1")]], m[[paste0(col, "_r2")]])
      row[[paste0("ac2_", col)]]      <- r$ac2
      row[[paste0("obsagree_", col)]] <- r$obs_agree
    }
    irr_score_rows[[p$group]] <- row
  }

  item_tab  <- bind_rows(irr_item_rows)
  score_tab <- bind_rows(irr_score_rows)

  cat("\nMean kappa (all items) [95% CI]:\n")
  print(round(naive_ci(item_tab$kappa), 3))
  cat("Mean observed agreement (all items):", round(mean(item_tab$obs_agree), 1), "%\n")

  cat("\nSensitivity analysis - manually-rated items only (excludes the 7",
      "auto-prefilled items) [95% CI]:\n")
  print(round(naive_ci(item_tab$kappa_manual), 3))
  cat("Mean observed agreement (manually-rated items only):",
      round(mean(item_tab$obs_agree_manual), 1), "%\n")

  cat("\nMean Gwet's AC2 per domain [95% CI], with mean observed agreement:\n")
  for (col in c("rob", "d1", "d2", "d3", "d4", "d5")) {
    vals <- score_tab[[paste0("ac2_", col)]]
    cat(" ", col, ": AC2 ="); print(round(naive_ci(vals), 3))
    cat("     observed agreement =",
        round(mean(score_tab[[paste0("obsagree_", col)]]), 1), "%\n")
  }
}, error = function(e) {
  cat("Objective 1 requires the metapsyTools and irrCAC packages",
      "(see manuscript references). Skipped:", conditionMessage(e), "\n")
})

# =============================================================================
# Load the shared instrument-level dataset for Objectives 2-5
# =============================================================================

dat <- read.csv(file.path(data_dir, "analysis_data.csv"), na.strings = "") %>%
  mutate(
    across(c(rob, rob_study_lvl, rob_d1, rob_d2, rob_d3, rob_d4, rob_d5),
           ~ tolower(trimws(.x))),
    control       = factor(control,       levels = c("other", "waitlist")),
    psy_type      = factor(psy_type,      levels = c("other", "CBT")),
    delivery      = factor(delivery,      levels = c("individual", "group", "guided self-help", "other")),
    diagnosis_cat = factor(diagnosis_cat, levels = c("clinical diagnosis", "cut-off/subclinical")),
    income_cat    = factor(income_cat,    levels = c("high", "low/middle")),
    log_n         = log(total_n)
  )

# One row per study; where ratings differ across rows (comparisons/instruments),
# take the most conservative rating, as in the Methods
most_conservative <- function(x) {
  f <- factor(tolower(trimws(x)), levels = rob_levels, ordered = TRUE)
  if (all(is.na(f))) NA_character_ else as.character(max(f, na.rm = TRUE))
}
study_dat <- dat %>%
  group_by(study) %>%
  mutate(across(c(rob_study_lvl, rob_d1, rob_d2, rob_d3, rob_d4, rob_d5), most_conservative)) %>%
  ungroup() %>%
  distinct(study, .keep_all = TRUE)

# Objective 5 uses effect-size-level control and psy_type (as in the original analysis)
dat_es <- dat %>%
  mutate(control  = factor(control_es,  levels = c("other", "waitlist")),
         psy_type = factor(psy_type_es, levels = c("other", "CBT")))

# =============================================================================
# OBJECTIVE 2 - Prevalence of risk of bias
# =============================================================================

cat("\n\n=== OBJECTIVE 2: Prevalence ===\n")

score_cols <- c("rob_d1", "rob_d2", "rob_d3", "rob_d4", "rob_d5", "rob_study_lvl")

# Study-level prevalence per domain uses the same "most conservative rating
# across instruments/comparisons" rule used to derive the published overall
# study-level RoB score (Methods, "Overall Risk of Bias Score"), applied
# independently to each domain. A simple one-row-per-study dedup is NOT
# equivalent: a handful of studies carry more than one rating per domain
# (most notably Domain 4, rated per instrument rather than per study - see
# Table 1), and picking an arbitrary row understates their true risk.
conservative_per_study <- function(data, col) {
  data %>%
    mutate(.f = factor(.data[[col]], levels = rob_levels, ordered = TRUE)) %>%
    filter(!is.na(.f)) %>%
    group_by(study) %>%
    summarise(category = factor(as.character(max(.f)), levels = rob_levels), .groups = "drop")
}

prevalence <- bind_rows(lapply(score_cols, function(col) {
  conservative_per_study(dat, col) %>%
    count(category, .drop = FALSE) %>%
    mutate(domain = domain_labels[col], pct = round(100 * n / sum(n), 1))
}))
print(as.data.frame(prevalence %>% select(domain, category, n, pct)), row.names = FALSE)

# Sensitivity: trials published 2001 onwards (revised CONSORT statement)
prevalence_2001 <- bind_rows(lapply(score_cols, function(col) {
  conservative_per_study(dat %>% filter(year >= 2001), col) %>%
    count(category, .drop = FALSE) %>%
    mutate(domain = domain_labels[col], pct = round(100 * n / sum(n), 1))
}))

# =============================================================================
# OBJECTIVE 3 - Trends in risk of bias over time
# =============================================================================

cat("\n\n=== OBJECTIVE 3: Trends over time ===\n")

run_poisson_trend <- function(data, col, label) {
  counts <- data %>%
    filter(!is.na(.data[[col]]), !is.na(year)) %>%
    mutate(category = factor(.data[[col]], levels = rob_levels)) %>%
    count(year, category, .drop = FALSE)
  counts$category <- relevel(counts$category, ref = "high risk")
  m <- glm(n ~ scale(year) * category, data = counts, family = poisson(link = "log"))
  cat("\n---", label, "---\n")
  print(round(summary(m)$coefficients, 4))
  m
}

trend_models <- Map(run_poisson_trend,
                     list(study_dat, study_dat, study_dat, study_dat, study_dat, study_dat),
                     score_cols, domain_labels[score_cols])

# Sensitivity analysis restricted to post-2001 trials was not performed for
# this objective (restricting the sample would remove the baseline period
# needed to assess a trend); documented in Supplement F.

# =============================================================================
# OBJECTIVE 4 - Trial characteristics associated with risk of bias
# =============================================================================

cat("\n\n=== OBJECTIVE 4: Trial characteristics (proportional odds models) ===\n")

fit_ordinal <- function(data, col, label) {
  d <- data %>% mutate(rob_ord = factor(.data[[col]], levels = rob_levels, ordered = TRUE))
  m <- polr(rob_ord ~ control + psy_type + delivery + diagnosis_cat + income_cat + log_n,
            data = d, Hess = TRUE)
  ctab <- coef(summary(m))
  pred_rows <- !grepl("\\|", rownames(ctab))
  p <- pnorm(abs(ctab[, "t value"]), lower.tail = FALSE) * 2
  out <- data.frame(
    domain = label, predictor = rownames(ctab)[pred_rows],
    OR    = round(exp(ctab[pred_rows, "Value"]), 2),
    CI_lo = round(exp(ctab[pred_rows, "Value"] - 1.96 * ctab[pred_rows, "Std. Error"]), 2),
    CI_hi = round(exp(ctab[pred_rows, "Value"] + 1.96 * ctab[pred_rows, "Std. Error"]), 2),
    p     = round(p[pred_rows], 4), n = nobs(m)
  )
  cat("\n---", label, "---\n"); print(out, row.names = FALSE)
  out
}

ordinal_results <- bind_rows(Map(fit_ordinal,
                                  list(study_dat, study_dat, study_dat, study_dat, study_dat, study_dat),
                                  score_cols, domain_labels[score_cols]))

# Sensitivity: trials published 2001 onwards
study_dat_2001 <- study_dat %>% filter(year >= 2001)
ordinal_results_2001 <- bind_rows(Map(fit_ordinal,
                                       list(study_dat_2001, study_dat_2001, study_dat_2001,
                                            study_dat_2001, study_dat_2001, study_dat_2001),
                                       score_cols, domain_labels[score_cols]))

# =============================================================================
# OBJECTIVE 5 - Association between risk of bias and effect size estimates
# =============================================================================

cat("\n\n=== OBJECTIVE 5: Risk of bias and effect size (CHE models) ===\n")

r_within <- 0.6
es_cols  <- c("rob", "rob_d1", "rob_d2", "rob_d3", "rob_d4", "rob_d5")

fit_che <- function(data, col, label, adjusted = FALSE) {
  d <- data %>%
    filter(!is.na(.data[[col]]), !is.na(.g_), !is.na(.g_se_)) %>%
    mutate(rob_fac = factor(.data[[col]], levels = rob_levels), vi = .g_se_^2)
  if (adjusted) d <- d %>% filter(!is.na(control))
  if (nrow(d) < 10) return(NULL)

  V   <- vcalc(vi = d$vi, cluster = d$study, obs = seq_len(nrow(d)), rho = r_within, nearpd = TRUE)
  fml <- if (adjusted) .g_ ~ rob_fac + control - 1 else .g_ ~ rob_fac - 1
  d   <- d %>% mutate(.row_id = seq_len(nrow(d)))

  m <- tryCatch(
    rma.mv(.g_, V, mods = fml, random = list(~ 1 | study / .row_id),
           data = d, method = "REML", test = "t", sparse = TRUE),
    error = function(e) { cat("  ERROR (", label, "):", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)

  cats     <- levels(d$rob_fac)[levels(d$rob_fac) %in% unique(d$rob_fac)]
  n_cats   <- length(cats)
  n_coef   <- length(coef(m))
  pooled   <- bind_rows(lapply(seq_len(n_cats), function(i) {
    nm   <- rep(0, n_coef); nm[i] <- 1
    pred <- predict(m, newmods = nm)
    data.frame(domain = label, adjusted = adjusted, category = cats[i],
               g = round(pred$pred, 3), ci_lo = round(pred$ci.lb, 3), ci_hi = round(pred$ci.ub, 3))
  }))

  m_rob <- tryCatch(robust(m, cluster = d$study, clubSandwich = TRUE), error = function(e) NULL)
  vcov_use <- if (!is.null(m_rob)) vcov(m_rob) else vcov(m)
  ddf_use  <- if (!is.null(m_rob)) min(m_rob$ddf) else min(m$ddf)
  coefs    <- coef(m)

  contrasts <- NULL
  if (n_cats >= 2) {
    idx_pairs <- combn(seq_len(n_cats), 2)
    contrasts <- bind_rows(lapply(seq_len(ncol(idx_pairs)), function(j) {
      i1 <- idx_pairs[1, j]; i2 <- idx_pairs[2, j]
      est <- coefs[i1] - coefs[i2]
      se  <- sqrt(vcov_use[i1, i1] + vcov_use[i2, i2] - 2 * vcov_use[i1, i2])
      data.frame(domain = label, adjusted = adjusted,
                 contrast = paste(cats[i1], "vs", cats[i2]),
                 delta_g = round(est, 3), se = round(se, 3),
                 ci_lo = round(est - qt(.975, ddf_use) * se, 3),
                 ci_hi = round(est + qt(.975, ddf_use) * se, 3))
    }))
  }
  list(pooled = pooled, contrasts = contrasts)
}

es_results <- list()
for (i in seq_along(es_cols)) {
  for (adj in c(FALSE, TRUE)) {
    res <- fit_che(dat_es, es_cols[i], domain_labels[gsub("^rob$", "rob_study_lvl", es_cols[i])], adjusted = adj)
    if (!is.null(res)) {
      cat("\n---", domain_labels[gsub("^rob$", "rob_study_lvl", es_cols[i])],
          ifelse(adj, "(adjusted)", "(unadjusted)"), "---\n")
      print(res$pooled, row.names = FALSE)
      print(res$contrasts, row.names = FALSE)
      es_results[[paste0(es_cols[i], if (adj) "_adj" else "_unadj")]] <- res
    }
  }
}

# Sensitivity: trials published 2001 onwards
dat_2001 <- dat_es %>% filter(year >= 2001)

es_results_2001 <- list()
for (i in seq_along(es_cols)) {
  for (adj in c(FALSE, TRUE)) {
    lab <- domain_labels[gsub("^rob$", "rob_study_lvl", es_cols[i])]
    res <- fit_che(dat_2001, es_cols[i], lab, adjusted = adj)
    if (!is.null(res)) {
      cat("\n--- 2001+:", lab, ifelse(adj, "(adjusted)", "(unadjusted)"), "---\n")
      print(res$pooled, row.names = FALSE)
      print(res$contrasts, row.names = FALSE)
      es_results_2001[[paste0(es_cols[i], if (adj) "_adj" else "_unadj")]] <- res
    }
  }
}

# Exploratory moderation: RoB dichotomized (low risk vs. higher risk) x
# study characteristics
cat("\n\n=== Exploratory moderation (RoB dichotomized x study characteristics) ===\n")
dat_mod <- dat_es %>%
  mutate(rob_dich = factor(ifelse(rob == "low risk", "low risk", "higher risk"),
                            levels = c("low risk", "higher risk")))

for (mod in c("control", "psy_type", "delivery", "diagnosis_cat", "income_cat", "log_n")) {
  d <- dat_mod %>% filter(!is.na(rob_dich), !is.na(.data[[mod]]), !is.na(.g_), !is.na(.g_se_)) %>%
    mutate(vi = .g_se_^2, .row_id = seq_len(n()))
  V   <- vcalc(vi = d$vi, cluster = d$study, obs = seq_len(nrow(d)), rho = r_within, nearpd = TRUE)
  fml <- as.formula(paste(".g_ ~ rob_dich *", mod))
  m   <- tryCatch(
    rma.mv(.g_, V, mods = fml, random = list(~ 1 | study / .row_id),
           data = d, method = "REML", test = "t", sparse = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) next
  ctab <- coef(summary(m))
  cat("\nModerator:", mod, "\n")
  print(round(ctab[grepl(":", rownames(ctab)), c("estimate", "se", "pval")], 4))
}

cat("\n\nDone.\n")

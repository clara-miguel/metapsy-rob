# Codebook

## nalysis_data.csv

One row per effect size (one depression instrument in one comparison): 514 trials, 941 effect sizes.

| Variable | Description | Values |
|---|---|---|
| `study` | Trial identifier (first author and year) | text |
| `year` | Publication year | integer |
| `rob` | Overall risk of bias score at the instrument level (used for Objective 5) | low risk / some concerns / high risk |
| `rob_study_lvl` | Overall risk of bias score at the study level, i.e., the most conservative rating across instruments (used for Objectives 2–4) | low risk / some concerns / high risk |
| `rob_d1` | Domain 1: randomization process | low risk / some concerns / high risk |
| `rob_d2` | Domain 2: deviations from intended interventions | same |
| `rob_d3` | Domain 3: missing outcome data | same |
| `rob_d4` | Domain 4: measurement of the outcome (rated per instrument) | same |
| `rob_d5` | Domain 5: selective outcome reporting | same |
| `.g_` | Hedges' g, psychotherapy vs. control at post-treatment (positive values favor psychotherapy) | numeric |
| `.g_se_` | Standard error of Hedges' g | numeric |
| `control_es` | Control condition of the comparison (effect-size level) | waitlist / other |
| `psy_type_es` | Psychotherapy type of the comparison (effect-size level) | CBT / other |
| `control` | Control condition at the study level (one value per trial, repeated across its rows) | waitlist (414 rows) / other (527 rows) |
| `psy_type` | Psychotherapy type at the study level (one value per trial, repeated across its rows) | CBT (512 rows) / other (429 rows) |
| `delivery` | Delivery format | individual (305) / group (249) / guided self-help (200) / other (187) |
| `diagnosis_cat` | Inclusion criterion for depression | clinical diagnosis (448) / cut-off/subclinical (493) |
| `income_cat` | World Bank income level of the country at the time of publication | high (772) / low/middle (166); missing for 3 rows (2 multinational trials) |
| `total_n` | Total sample size across arms | integer, 14–1,899, no missing values |

Notes:
- Domain ratings are used at the study level after taking the most conservative rating across a study's rows (done in the analysis script).
- Domains 1, 2, 3 and 5 are rated at the study or comparison level and Domain 4 per instrument, so ratings can differ between rows of the same trial.
- Overall rule: all five domains low risk = low risk; any domain high risk, or three or more domains with some concerns = high risk; otherwise some concerns.

## shared_data/irr_ratings.csv

Independent item-level ratings for the inter-rater reliability subsample (237 trials, 10 rater pairs, 16 raters). One row per rater per trial.

| Variable | Description |
|---|---|
| `group` | Rater pair |
| `rater` | Rater initials |
| `study` | Trial identifier (same as in `analysis_data.csv`) |
| `d1_1` to `d1_4` | Items 1–4 (Domain 1) |
| `d2_5` to `d2_9` | Items 5–9 (Domain 2) |
| `d3_10` to `d3_14` | Items 10–14 (Domain 3) |
| `d4_15` to `d4_18` | Items 15–18 (Domain 4) |
| `d5_19`, `d5_21` to `d5_24` | Items 19 and 21–24 (Domain 5) |

Response coding: Yes/PY (yes or probably yes), No/PN (no or probably no), NI (no information); blank = no rating.

Notes:
- Seven items are pre-filled from extracted data and are identical for both raters of a pair: d1_3, d1_4, d3_10, d3_12, d3_13, d4_16 and d4_17.
- The full questionnaire is in the supplements of the manuscript.

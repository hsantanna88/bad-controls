#' ============================================================================
#' Difference-in-Differences with Bad Controls: A Practical Guide
#' ============================================================================
#'
#' Authors: Carolina Caetano, Brantly Callaway, Stroud Payne, Hugo Sant'Anna
#'
#' This guide demonstrates:
#'   1. What the "bad controls" problem is in DiD
#'   2. How standard approaches fail
#'   3. How our proposed methods solve it
#'
#' We use a simulated dataset where the true ATT is known (ATT = 1.0),
#' so you can verify that the estimators recover the correct answer.
#' ============================================================================

library(badcontrols)
library(pte)
library(ggplot2)
library(dplyr)

set.seed(20240301)

# =============================================================================
# PART 1: THE BAD CONTROLS PROBLEM
# =============================================================================
#
# In difference-in-differences, researchers often want to include covariates
# in the parallel trends assumption:
#
#   E[Y_t(0) - Y_{t-1}(0) | X_t, D=1] = E[Y_t(0) - Y_{t-1}(0) | X_t, D=0]
#
# The problem: if X_t is measured AFTER treatment, then X_t is itself affected
# by treatment. Conditioning on it introduces post-treatment selection bias.
#
# Three things can go wrong:
#   (a) Including X_t directly introduces bias (bad control)
#   (b) Dropping X_t entirely may violate parallel trends
#   (c) Using only X_{t-1} may not fully solve the problem
# =============================================================================


# =============================================================================
# PART 2: SIMULATE DATA
# =============================================================================

sim <- simulate_bad_controls(
  n = 2000,
  T_max = 4,
  direct_att = 0.2,
  x_effect_on_y = 1.0,
  treatment_effect_on_x = 0.8
)

sim_data <- sim$data
true_att <- sim$true_att

cat("========================================\n")
cat("True TOTAL ATT =", true_att, "\n")
cat("  = direct_att (0.2) + beta * theta (1.0 * 0.8 = 0.8)\n")
cat("========================================\n\n")


# =============================================================================
# PART 3: VISUALIZE THE BAD CONTROLS PROBLEM
# =============================================================================

x_means <- sim_data %>%
  mutate(treated_ever = ifelse(G > 0, "Treated", "Never-Treated")) %>%
  group_by(treated_ever, period) %>%
  summarize(mean_X = mean(X), .groups = "drop")

p1 <- ggplot(x_means, aes(x = period, y = mean_X, color = treated_ever)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 2.5, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "The Bad Control: X diverges after treatment",
    subtitle = "Treatment causally shifts X, making post-treatment X endogenous",
    x = "Period", y = "Mean of X", color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p1)


# =============================================================================
# PART 4: COMPARE ESTIMATORS
# =============================================================================

# Method 0: Naive DiD (no covariates) -- biased
cat("\n--- Method 0: Naive DiD (no covariates) ---\n")
res_naive <- pte::pte_default(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim_data, d_outcome = TRUE, est_method = "reg"
)
att0 <- extract_att(res_naive)
cat("ATT:", round(att0$att, 4), "(SE:", round(att0$se, 4), ")\n")

# Method 1: Include X_t directly (bad control) -- biased
cat("\n--- Method 1: Include X_t (BAD CONTROL) ---\n")
res_bad <- pte::pte_default(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim_data, d_outcome = TRUE,
  d_covs_formula = ~X, est_method = "reg"
)
att1 <- extract_att(res_bad)
cat("ATT:", round(att1$att, 4), "(SE:", round(att1$se, 4), ")\n")

# Method 2: Pre-treatment X only
cat("\n--- Method 2: Pre-treatment X only ---\n")
res_pretreat <- pte::pte_default(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim_data, d_outcome = TRUE,
  xformla = ~Z, est_method = "reg"
)
att2 <- extract_att(res_pretreat)
cat("ATT:", round(att2$att, 4), "(SE:", round(att2$se, 4), ")\n")

# Method 3: Imputation (our proposal)
cat("\n--- Method 3: Imputation ---\n")
res_impute <- bc_att_gt(
  yname = "Y", gname = "G", tname = "period", idname = "id",
  data = sim_data,
  bad_control_formula = ~X,
  xformla = ~Z,
  est_method = "imputation"
)
att3 <- extract_att(res_impute)
cat("ATT:", round(att3$att, 4), "(SE:", round(att3$se, 4), ")\n")

# Method 4: Doubly Robust ML (our proposal)
cat("\n--- Method 4: Doubly Robust ML ---\n")
if (requireNamespace("grf", quietly = TRUE)) {
  res_drml <- bc_att_gt(
    yname = "Y", gname = "G", tname = "period", idname = "id",
    data = sim_data,
    bad_control_formula = ~X,
    xformla = ~Z,
    est_method = "dr_ml"
  )
  att4 <- extract_att(res_drml)
  cat("ATT:", round(att4$att, 4), "(SE:", round(att4$se, 4), ")\n")
} else {
  cat("Skipped: install 'grf' package.\n")
  att4 <- list(att = NA, se = NA)
}


# =============================================================================
# PART 5: COMPARISON TABLE
# =============================================================================

cat("\n=============================================================\n")
cat("COMPARISON OF ESTIMATORS (True ATT =", true_att, ")\n")
cat("=============================================================\n")

results <- data.frame(
  Method = c(
    "0: Naive DiD",
    "1: Bad Control",
    "2: Pre-treatment X",
    "3: Imputation",
    "4: DR/ML"
  ),
  ATT = c(att0$att, att1$att, att2$att, att3$att, att4$att),
  SE = c(att0$se, att1$se, att2$se, att3$se, att4$se)
)
results$Bias <- results$ATT - true_att

print(results[, c("Method", "ATT", "Bias", "SE")], row.names = FALSE)

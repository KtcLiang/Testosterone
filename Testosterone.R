# Meta-Analysis: PFSF Sexual Desire Domain
# Testosterone Therapy in Surgically Menopausal Women

# Load required packages
install.packages("metafor")
install.packages("ggplot2")
library(metafor)
library(ggplot2)

# ---- 1. Input Data ----
lsmean_diff <- c(5.3, 6.28, 5.12, 10.45)
se_diff <- c(2.70, 1.90, 1.49, 4.49)
ci_lower <- round(lsmean_diff-se_diff*qnorm(.975),2)
ci_upper <-round(lsmean_diff+se_diff*qnorm(.975),2)
meta_data <- data.frame(
  study = c("Braunstein 2005", "Buster 2005", "Simon 2005", "Davis 2006"),
  
  # Sample sizes
  n_trt = c(110, 267, 283, 37),
  n_ctrl = c(119, 266, 279, 40),
  
  # LS mean changes from baseline (ANCOVA-adjusted)
  lsmean_trt = c(13.7, 10.57, 11.06, 16.43),
  lsmean_ctrl = c(8.4, 4.29, 5.94, 5.98),
  
  # LS mean difference (treatment - control)
  lsmean_diff = lsmean_diff,
  
  # Standard errors
  se_diff = se_diff,
  
  # Two-sided P-values
  pval = c(0.05, 0.001, 0.0006, 0.02),
  # 95% Confidence intervals
  ci_lower=ci_lower,
  ci_upper=ci_upper,
  
  # Baseline PFSF Sexual Desire (pooled across arms)
  baseline_trt = c(20.9, 21.6, 19.8, 23.1),
  baseline_ctrl = c(20.9, 21.6, 20.8, 19.9),
  baseline_pool = c(20.9, 21.6, 20.3, 21.5),
  # Age (study-level mean, years)
  age_mean_trt = c(49.6, 48.3, 49.2, 51),
  age_mean_ctrl = c(48.5, 49.5, 48.9, 49.3),
  age_mean_pool = c(49.0, 48.9, 49.1, 50.1),
  # Time in relationship (years)
  time_relation_trt = c(17.2, 17.5, 19.7, 21.1),
  time_relation_ctrl = c(18.1, 18.9, 18.6, 21.1),
  time_relation_pool = c(17.7, 18.2, 19.2, 21.1),
  # Estrogen route
  estrogen_route = c("oral", "most oral", "most oral", "transdermal")
)

# Calculate variance
meta_data$var_diff <- meta_data$se_diff^2

# Display raw data
print("=== RAW DATA ===")
print(meta_data)


# ---- 2. Descriptive Summary ----
print("\n=== DESCRIPTIVE SUMMARY ===")
summary_stats <- data.frame(
  Study = meta_data$study,
  N_Total = meta_data$n_trt + meta_data$n_ctrl,
  LS_Mean_Diff = meta_data$lsmean_diff,
  SE = round(meta_data$se_diff, 2),
  CI_95 = paste0("(", meta_data$ci_lower, ", ", meta_data$ci_upper, ")"),
  P_value = meta_data$pval,
  Baseline = meta_data$baseline_pool,
  Age = meta_data$age_mean_pool,
  Time_in_Relation = meta_data$time_relation_pool,
  estrogen_route=meta_data$estrogen_route
)
print(summary_stats)

# ---- 3. Random-Effects Meta-Analysis (Primary) ----
print("\n=== PRIMARY RANDOM-EFFECTS META-ANALYSIS ===")
# re_model <- rma(yi = lsmean_diff,
#                 sei = se_diff,
#                 data = meta_data,
#                 method = "DL",  # DerSimonian-Laird estimator
#                 slab = study)
re_model <- rma(yi = lsmean_diff,
                sei = se_diff,
                data = meta_data,
                method = "REML",  # RMEL estimator for heterogeneity
                slab = study)
print(re_model)

# Extract key results
cat("\n--- Pooled Estimate ---\n")
cat("LS Mean Difference: ", round(re_model$beta[1], 2), "\n")
cat("95% CI: (", round(re_model$ci.lb, 2), ", ", round(re_model$ci.ub, 2), ")\n")
cat("P-value: ", format.pval(re_model$pval, eps = 0.001), "\n")

#names(re_model)
cat("\n--- Heterogeneity ---\n")
cat("Tau²: ", round(re_model$tau2, 2), "\n")
cat("I²: ", round(re_model$I2, 1), "%\n")
Q_df   <- re_model$k - re_model$p
cat("Q-statistic: ", round(re_model$QE, 2), " (df=", Q_df, ", P=",
    format.pval(re_model$QEp, eps = 0.001), ")\n")

# Prediction interval (For a random-effects model, predict(object) returns the
#estimated (average) outcome in the hypothetical population of studies from which
#the set of studies included in the meta-analysis are assumed to be a random selection.
#This is the same as the estimated intercept in the random-effects model (i.e. estimated mu hat)
pred_int <- predict(re_model)
cat("\n--- 95% Prediction Interval ---\n")
cat("Future study range: (", round(pred_int$pi.lb, 2), ", ",
    round(pred_int$pi.ub, 2), ")\n")

# ---- 4. Meta-Regression with Covariates as moderator ----
print("\n=== META-REGRESSION WITH COVARIATES ===")
mr_model1 <- rma(yi = lsmean_diff,
                 sei = se_diff,
                 mods = ~ baseline_pool,
                 #mods = ~ age_mean_pool,
                 #mods = ~ time_relation_pool,
                 #mods = ~ estrogen_route,
                 data = meta_data,
                 method = "REML")  # Restricted maximum likelihood
mr_model2 <- rma(yi = lsmean_diff,
                 sei = se_diff,
                 #mods = ~ baseline_pool,
                 mods = ~ age_mean_pool,
                 ##mods = ~ time_relation_pool,
                 #mods = ~ estrogen_route,
                 data = meta_data,
                 method = "REML")
mr_model3 <- rma(yi = lsmean_diff,
                 sei = se_diff,
                 #mods = ~ baseline_pool,
                 #mods = ~ age_mean_pool,
                 mods = ~ time_relation_pool,
                 #mods = ~ estrogen_route,
                 data = meta_data,
                 method = "REML")
mr_model <- rma(yi = lsmean_diff,
                sei = se_diff,
                #mods = ~ baseline_pool,
                #mods = ~ age_mean_pool,
                #mods = ~ time_relation_pool,
                mods = ~ estrogen_route,
                data = meta_data,
                method = "REML")
print(mr_model1)
print(mr_model2)
print(mr_model3)
print(mr_model)
#For all four mr_model results--only for mr_model, intercept is significant greater than 0 when regress on estrogen_route,
#which means y is not 0 regardless which estrogen_route used.

# Interpret moderator effects
cat("\n--- Moderator Interpretation 1---\n")
coef_table <- data.frame(
  Moderator = c("Baseline PFSF Desire"),
  Coefficient = round(coef(mr_model1)[-1], 3),
  SE = round(mr_model1$se[-1], 3),
  P_value = format.pval(mr_model1$pval[-1], eps = 0.001)
)
print(coef_table)
cat("\nR² (variance explained): ", round(mr_model1$R2, 1), "%\n")

cat("\n--- Moderator Interpretation 2---\n")
coef_table <- data.frame(
  Moderator = c("Baseline Age"),
  Coefficient = round(coef(mr_model2)[-1], 3),
  SE = round(mr_model2$se[-1], 3),
  P_value = format.pval(mr_model2$pval[-1], eps = 0.001)
)
print(coef_table)
cat("\nR² (variance explained): ", round(mr_model2$R2, 1), "%\n")


cat("\n--- Moderator Interpretation 3---\n")
coef_table <- data.frame(
  Moderator = c("Baseline Time in Relationship"),
  Coefficient = round(coef(mr_model3)[-1], 3),
  SE = round(mr_model3$se[-1], 3),
  P_value = format.pval(mr_model3$pval[-1], eps = 0.001)
)
print(coef_table)
cat("\nR² (variance explained): ", round(mr_model3$R2, 1), "%\n")

cat("\n--- Moderator Interpretation 4---\n")
coef_table <- data.frame(
  Moderator = c("Estrogen Route"),
  Coefficient = round(coef(mr_model)[-1], 3),
  SE = round(mr_model$se[-1], 3),
  P_value = format.pval(mr_model$pval[-1], eps = 0.001)
)
print(coef_table)
cat("\nR² (variance explained): ", round(mr_model$R2, 1), "%\n")

# ---- 5. Forest Plot ----
print("\n=== GENERATING FOREST PLOT ===")

#weight <-weights(re_model)
# Save as PDF (optional)
pdf("forest_plot_pfsf_desire.pdf", width = 12, height = 12)

forest(re_model,showweights=F,
       xlim = c(-8, 25),
       at = c(0, 5, 10, 15, 20),
       xlab = "LS Mean Difference in PFSF Sexual Desire (0-100 scale)",
       slab = paste(meta_data$study, " (N=", meta_data$n_trt + meta_data$n_ctrl, ")",  "\n ", round(meta_data$baseline_pool,1), ", ", round(meta_data$age_mean_pool,1), ", ",
                    round(meta_data$time_relation_pool,1), ", ", meta_data$estrogen_route, sep=""),
       #"\n(", "N=",
       #meta_data$n_trt + meta_data$n_ctrl, ")", sep = ""),
       #ilab = cbind(meta_data$baseline_pool, meta_data$age_mean_pool),
       #ilab=c("Weight", "Eestimate (Prediction CIs)"),
       #ilab.xpos = c(-6, -8),
       #ilab.xpos = c(15, 15),
       cex = 1,
       header = c("Study (Total N): Testosterone 300 µg/d vs. Pbo", "Est (Pred CIs)"),
       mlab = "Pooled Estimate (Random-Effects)",
       col = "darkblue",
       border = "darkblue",
       refline = 0,
       addpred = TRUE,  # Add prediction interval
       lwd = 2)
#par("usr") get the limit for x, y
# Add column headers
#text(c(20, 15), re_model$k + 4, c("Baseline", "Age", "Time in Relation", "E Route"), cex = 0.9, font = 2)
#text(c(0, 5.5), pos=4, c("Baseline PFSF, Age, Time_in_Relation, Route"), cex = 0.7, font = 2)
#text(c(-8, 6.5), pos=4, c("Baseline PFSF, Age, Time_in_Relation, Route"), cex = 0.7, font = 2)
#text(c(0, 7.5), pos=4, c("Baseline PFSF, Age, Time_in_Relation, Route"), cex = 0.7, font = 2)
text(-8, re_model$k+1.5, pos=4, c("Baseline PFSF, Age, Time_in_Relation, Route"), cex = 0.9, font = 2)

# Add reference line and labels
abline(v = 0, lty = 2, col = "gray50")
text(1, -1.5, pos = 4, cex = 0.8, font = 3,
     "← Favors Pbo          Favors Testosterone →")
text(0, -2.5, pos = 4, cex = 0.75, font = 3,
     paste0("Heterogeneity: I² = ", round(re_model$I2, 1), "%, τ² = ",
            round(re_model$tau2, 2), ", Q = ", round(re_model$QE, 2),
            " (P = ", format.pval(re_model$QEp, eps = 0.01), ")"))

# Close PDF device if opened
dev.off()

print("Forest plot generated successfully.")

# ---- 6. Sensitivity Analyses ----
print("\n=== SENSITIVITY ANALYSES ===")

# 6a. Leave-one-out analysis
print("\n--- Leave-One-Out Analysis ---")
loo_results <- leave1out(re_model)
print(loo_results)
#all leave-one out models still has significant treatment effect

# 6b. Influence diagnostics
print("\n--- Influential Studies ---")
inf <- influence(re_model)
print(inf)
#Sensitivity analysis to assess influence of each individual study (completes the leave-one-out assessment)
#rstudent         dffits cook.d  cov.r   tau2.del QE.del    hat  weight    dfbs  inf
#Braunstein 2005  -0.1959 -0.0823 0.0068 1.1765   0.0000 1.3404 0.1500 15.0046 -0.0823
#Buster 2005       0.3104  0.2046 0.0419 1.4347   0.0000 1.2825 0.3030 30.3001  0.2046
#Simon 2005       -0.6291 -0.6200 0.3844 1.9712   0.0000 0.9830 0.4927 49.2696 -0.6200
#Davis 2006        1.0677  0.2557 0.0654 1.0574   0.0000 0.2387 0.0543  5.4257  0.2557
#Braunstein 2005 and Buster 2005 not influential, Simon 2005 has the most leverage (49.3% weight)
#due to its precision (small SE) and it drives the pooled estimate the most.
#Davis 2006 is likely outlier, the primary contributor to heterogeneity (tau2.del = 0 when removed).

# ---- 7. Publication Bias Assessment ----
print("\n=== PUBLICATION BIAS ASSESSMENT ===")

# Funnel plot
print("Generating funnel plot...")
funnel(re_model,
       xlab = "LS Mean Difference",
       main = "Funnel Plot: PFSF Sexual Desire")

# Egger's test (limited power with k=4)
print("\n--- Egger's Regression Test ---")
egger_test <- regtest(re_model, model = "lm")
print(egger_test)

#Symmetry Assessment
#The funnel appears slightly asymmetric — Davis 2006 sits at the bottom-right (large effect, large SE), with no corresponding study at the bottom-left (large SE, small/negative effect).
#However, with only k = 4 studies, it is very difficult to draw meaningful conclusions about symmetry. Funnel plots are generally unreliable with fewer than 10 studies.
#The Inverted Funnel Shape
#The two most precise studies (Simon 2005, Buster 2005) cluster tightly near the top around the pooled estimate (~5.79), forming the narrow tip of the funnel.
#The less precise studies (Braunstein 2005, Davis 2006) spread out toward the bottom — this is expected behavior.
#sample size too small to draw conclusion that there is publication bias. The assymetry caused by Davis study likely not due to publication bias.
#Davis 2006 is an outlier to the right             Could suggest publication bias (missing small studies with null/negative results on the left) OR could reflect genuine heterogeneity (Davis used transdermal estrogen, unlike the others)
#No studies in the bottom-left quadrant     Potentially concerning for publication bias, but sample is too small to conclude
#Top studies cluster symmetrically Reassuring — the most precise studies agree well
#Egger's regression test examines whether there is a systematic relationship between effect size and study precision (standard error).
#Although test result not significant, the test is unlikely to detect it with so few studies even if publication bias exists, .

# ---- 8. Export Results ----
print("\n=== EXPORTING RESULTS ===")

# Export summary table
write.csv(summary_stats, "meta_analysis_summary.csv", row.names = FALSE)
cat("Summary table saved to: meta_analysis_summary.csv\n")

# Export model results
sink("meta_analysis_results.txt")
cat("=== PRIMARY RANDOM-EFFECTS MODEL ===\n")
print(re_model)
cat("\n=== META-REGRESSION MODEL ===\n")
print(mr_model)
sink()
cat("Model results saved to: meta_analysis_results.txt\n")

print("\n=== ANALYSIS COMPLETE ===")
######## TMT sRBD hyp vs controls longitudinal variance analysis ########
library(markdown)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(e1071)
library(gridExtra)
library(broom.mixed)
library(boot)
library(readr)
library(ggdist) 
library(nlme)
library(lme4)
library(lmerTest)
library(nlme) 
library(tidybayes)
library(moments)


matched_sRBD_hyp_control_TMT_data <- read_csv("TMT/sRBD hyp vs controls/matched_sRBD_hyp_control_TMT_data.csv")
View(matched_sRBD_hyp_control_TMT_data)


table(matched_sRBD_hyp_control_TMT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
matched_sRBD_hyp_control_TMT_data <- matched_sRBD_hyp_control_TMT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))

matched_sRBD_hyp_control_TMT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)

matched_sRBD_hyp_control_TMT_data %>%
  filter(Code == "Control", PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline != "Unknown") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)


# Define the columns of interest
tmt_cols <- c(
  "tmt_trailAminustrailB_exefcn_final_4", "tmt_trailAminustrailB_exefcn_final_5", "tmt_trailAminustrailB_exefcn_final_7",
  "tmt_trailAminustrailB_exefcn_final_8", "tmt_trailB_reactiontime_4", "tmt_trailB_reactiontime_5", "tmt_trailB_reactiontime_7",
  "tmt_trailB_reactiontime_8", "tmt_trailA_reactiontime_4", "tmt_trailA_reactiontime_5", "tmt_trailAminustrailB_exefcn_final_7",
  "tmt_trailA_reactiontime_8"
)

# Remove participants with all test data missing
matched_sRBD_hyp_control_TMT_data2 <- matched_sRBD_hyp_control_TMT_data %>%
  filter(
    rowSums(is.na(select(., all_of(tmt_cols)))) <
      ncol(select(., all_of(tmt_cols)))
  )


write.csv(matched_sRBD_hyp_control_TMT_data2, "TMT/sRBD hyp vs controls/matched_sRBD_hyp_control_TMT_data2.csv",
          row.names = FALSE)



######## Exploratory data analysis #########
# List of your TMT columns
cols_tmt <- c(
  "tmt_trailAminustrailB_exefcn_final_4", "tmt_trailAminustrailB_exefcn_final_5", "tmt_trailAminustrailB_exefcn_final_7",
  "tmt_trailAminustrailB_exefcn_final_8", "tmt_trailB_reactiontime_4", "tmt_trailB_reactiontime_5", "tmt_trailB_reactiontime_7",
  "tmt_trailB_reactiontime_8", "tmt_trailA_reactiontime_4", "tmt_trailA_reactiontime_5", "tmt_trailAminustrailB_exefcn_final_7",
  "tmt_trailA_reactiontime_8")

#check columns
names(matched_sRBD_hyp_control_TMT_data2)[!sapply(matched_sRBD_hyp_control_TMT_data2, is.numeric)]


#make columns numeric
cols_to_convert <- c(
  "tmt_trailAminustrailB_exefcn_final_4", "tmt_trailAminustrailB_exefcn_final_5",
  "tmt_trailAminustrailB_exefcn_final_7", "tmt_trailAminustrailB_exefcn_final_8",
  "tmt_trailB_reactiontime_4", "tmt_trailB_reactiontime_5",
  "tmt_trailB_reactiontime_7", "tmt_trailB_reactiontime_8",
  "tmt_trailA_reactiontime_4", "tmt_trailA_reactiontime_5",
  "tmt_trailA_reactiontime_7",
  "tmt_trailA_reactiontime_8"
)

matched_sRBD_hyp_control_TMT_data2[cols_to_convert] <- lapply(matched_sRBD_hyp_control_TMT_data2[cols_to_convert], as.numeric)

#remove participants with ALL test data missing 
matched_sRBD_hyp_control_TMT_data2 <- matched_sRBD_hyp_control_TMT_data2 %>%
  filter(
    rowSums(is.na(select(., all_of(cols_to_convert)))) <
      ncol(select(., all_of(cols_to_convert)))
  )


# Pivot to long format but keep all non test columns
sRBD_hyp_TMT_long <- matched_sRBD_hyp_control_TMT_data2 %>%
  pivot_longer(
    cols = all_of(cols_tmt),                 
    names_to = c("Measure", "Session"),
    names_pattern = "(.*)_(\\d+)$",          
    values_to = "Value"
  ) %>%
  mutate(Session = as.numeric(Session)) %>%
  group_by(PID2, Measure) %>%
  arrange(PID2, Measure, Session) %>%
  mutate(RelativeVisit = row_number()) %>%   # integer per participant
  ungroup()


# check how many timepoints
sRBD_hyp_TMT_long %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)

sRBD_hyp_TMT_long %>%
  filter(!is.na(Value)) %>%                       
  group_by(Measure, PID2) %>%
  summarise(n_obs = n(), .groups = "drop") %>%   
  filter(n_obs == 4) %>%                       
  count(Measure, name = "n_IDs_complete_4tp")


#Variability-eligible dataset (maximising sample)
sRBD_hyp_TMT_long_all <- sRBD_hyp_TMT_long %>%
  filter(!is.na(Value)) %>%
  group_by(PID2, Code, Measure) %>%
  filter(n() >= 2) %>%                 
  arrange(Session, .by_group = TRUE) %>%
  mutate(RelativeVisit = row_number()) %>%  
  ungroup()

sRBD_hyp_TMT_long_all %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)


#check if RT is actually 0
sRBD_hyp_TMT_long_all %>%
  filter(
    Measure %in% c("tmt_trailA_reactiontime", "tmt_trailB_reactiontime"),
    Value == 0
  ) %>%
  select(PID2, Code, Session, RelativeVisit, Measure, Value) %>%
  arrange(PID2, Session)
# no zeros

# Horizontal raincloud plots per visit
visits <- sort(unique(sRBD_hyp_TMT_long_all$RelativeVisit))

for (v in visits) {
  
  df_v <- sRBD_hyp_TMT_long_all %>%
    filter(RelativeVisit == v) %>%
    filter(is.finite(Value))   # keep only finite values (prevents blank panels)
  
  if (nrow(df_v) == 0) next
  
  p <- df_v %>%
    ggplot(aes(y = Code, x = Value, fill = Code, color = Code)) +
    stat_halfeye(
      adjust = 0.5, width = 0.6, .width = 0,
      justification = -0.3, orientation = "horizontal"
    ) +
    geom_boxplot(
      width = 0.2, outlier.shape = NA,
      alpha = 0.5, position = position_nudge(y = 0.3)
    ) +
    geom_jitter(height = 0.1, size = 1.5, alpha = 0.8) +
    facet_wrap(~Measure, scales = "free_x") +
    scale_fill_manual(values = c("Control" = "#1f3552",
                                 "sRBD_hyposmia" = "#69b3a2")) +
    scale_color_manual(values = c("Control" = "#1f3552",
                                  "sRBD_hyposmia" = "#69b3a2")) +
    labs(
      y = "Group",
      x = "Score",
      title = paste("Raincloud Plot of TMT Measures – Visit", v)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      plot.margin = margin(10, 30, 10, 10)  # <- adds space on the right
    ) +
    coord_cartesian(clip = "off")  # <- prevents cutting off last tick label
  
  print(p)
}


#check for outliers
sRBD_hyp_TMT_flagged <- sRBD_hyp_TMT_long_all %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
  ) %>%
  ungroup()


# Flag outliers and create summary table in one go
outlier_counts <- sRBD_hyp_TMT_long_all %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
  ) %>%
  ungroup() %>%
  group_by(Code) %>%
  summarise(
    total_outliers = sum(is_outlier, na.rm = TRUE),
    total_points = n(),
    perc_outliers = round(100 * total_outliers / total_points, 2)
  )

outlier_counts


# Identify outliers by Code, measure, and year
outlier_data <- sRBD_hyp_TMT_long_all %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
  ) %>%
  ungroup()


# Boxplot highlighting outliers
ggplot(outlier_data, aes(x = factor(RelativeVisit), y = Value, fill = Code)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.color = "red", alpha = 0.7) +
  facet_wrap(~Measure, scales = "free_y") +
  labs(
    title = "Distribution of TMT Scores by Year and Group",
    x = "Year",
    y = "Score Value",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


sRBD_hyp_TMT_long_all <- sRBD_hyp_TMT_long_all %>%
  mutate(RelativeVisit = as.numeric(RelativeVisit))


#check if data is normally distributed 
ggplot(sRBD_hyp_TMT_long_all %>% filter(Measure == "tmt_trailA_reactiontime"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Trail A Reaction Time by Group",
    x = "Trail A Reaction Time",
    y = "Density",
    fill = "Group"
  ) +
  xlim(0,3000) +
  theme_minimal(base_size = 14)


ggplot(sRBD_hyp_TMT_long_all %>% filter(Measure == "tmt_trailB_reactiontime"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Trail B Reaction Time by Group",
    x = "Trail B Reaction Time",
    y = "Density",
    fill = "Group"
  ) +
  xlim(0,6500) +
  theme_minimal(base_size = 14)


ggplot(sRBD_hyp_TMT_long_all %>% filter(Measure == "tmt_trailAminustrailB_exefcn_final"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Executive Function by Group",
    x = "Executive Function",
    y = "Density",
    fill = "Group"
  ) +
  xlim(-1050,3500) +
  theme_minimal(base_size = 14)


#test for normality
normality_results <- sRBD_hyp_TMT_long_all %>%
  group_by(Measure, PID2) %>%       # assume ID is your participant identifier
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%  # one value per person
  group_by(Measure) %>%
  summarise(
    n = sum(!is.na(Value)),
    p_value = tryCatch(
      if (n >= 3) shapiro.test(Value)$p.value else NA_real_,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(p_value) ~ "Not tested (too few values)",
      p_value > 0.05 ~ "Likely normal",
      TRUE ~ "Non-normal"
    )
  )

normality_results


#inspect visually
# Loop through each measure
for(m in unique(sRBD_hyp_TMT_long_all$Measure)) {
  data_m <- sRBD_hyp_TMT_long_all %>% filter(Measure == m, !is.na(Value))
  
  # Print skewness
  cat(m, "skewness:", round(skewness(data_m$Value), 2), "\n")
  
  # Plot histogram + QQ plot side by side
  gridExtra::grid.arrange(
    ggplot(data_m, aes(x = Value)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "#69b3a2", alpha = 0.6) +
      geom_density(color = "#1f3552", size = 1) +
      labs(title = paste(m, "- Histogram & Density"), x = "Value", y = "Density") +
      theme_minimal(base_size = 12),
    
    ggplot(data_m, aes(sample = Value)) +
      stat_qq() + stat_qq_line(color = "red") +
      labs(title = paste(m, "- QQ Plot"), x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal(base_size = 12),
    
    ncol = 2
  )
}


######### Variance Analysis ######### 
#Step 1a — Compute mean, SD, and sample size per group × year × measure
variance_summary <- sRBD_hyp_TMT_long_all %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value)),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(variance_summary)

#Step 1b — Bootstrap SD for confidence intervals (small sample size)
# Function to compute SD for bootstrapping
sd_fun <- function(data, indices) {
  return(sd(data[indices], na.rm = TRUE))
}

# Bootstrap SD function per group
bootstrap_sd_ci <- function(x, R = 1000, conf = 0.95) {
  boot_obj <- boot(data = x, statistic = sd_fun, R = R)
  ci <- boot.ci(boot_obj, type = "perc")$percent[4:5]  # lower and upper CI
  return(data.frame(sd = sd(x, na.rm = TRUE),
                    sd_lower = ci[1],
                    sd_upper = ci[2]))
}


# Apply to all groups
variance_boot <- sRBD_hyp_TMT_long_all %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value)),
    boot_results = list(bootstrap_sd_ci(Value)),
    .groups = "drop"
  ) %>%
  unnest(cols = c(boot_results))


print(variance_boot, n = Inf, width = Inf)

write.csv(variance_boot, "TMT/sRBD hyp vs controls/variance_boot_sRBD_hyp.csv", row.names = FALSE)


#Step 1c — Plot mean ± SD with bootstrap CI
# Combine mean and bootstrap SD
plot_data <- variance_summary %>%
  left_join(variance_boot, by = c("Measure", "RelativeVisit", "Code"))


# Get unique measures
measures <- unique(plot_data$Measure)

# Named vector for nice measure labels
measure_labels <- c(
  "tmt_trailAminustrailB_exefcn_final" = "Executive Function",
  "tmt_trailA_reactiontime" = "Trail A Reaction Time",
  "tmt_trailB_reactiontime" = "Trail B Reaction Time"
)


plot_data <- plot_data %>%
  mutate(Code = ifelse(Code == "sRBD_hyposmia", "pRBD with hyposmia", Code))


# Loop through measures and plot with subtitle
for (m in measures) {
  
  p <- plot_data %>%
    filter(Measure == m) %>%
    ggplot(aes(x = RelativeVisit, y = mean_value, color = Code, group = Code, fill = Code)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_value - sd_lower, ymax = mean_value + sd_upper),
                alpha = 0.20, color = NA) +
    labs(
      title = paste0("TMT – ", measure_labels[m], " Over Time (Mean ± SD)"),
      x = "Year",
      y = "Score",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic"),
      legend.position = "top"
    )
  
  print(p)
}


# Save the combined mean + SD + bootstrap CI
write.csv(plot_data, "TMT/sRBD hyp vs controls/TMT_mean_SD_bootstrap_sRBD_hyp.csv", row.names = FALSE)



#Step 2: subject-level variance (“wobblers”) + wilcoxon rank
calc_subject_variability <- function(long_df,
                                     measures,
                                     id_col   = "PID2",
                                     group_col = "Code",
                                     time_col = "RelativeVisit",
                                     value_col = "Value",
                                     exact_two_sided = FALSE) {
  
  stopifnot(all(c("Measure", id_col, group_col, time_col, value_col) %in% names(long_df)))
  stopifnot(all(measures %in% unique(long_df$Measure)))
  
  # Helper: compute SD, Range, MSSD for one subject (already ordered by time)
  subject_metrics <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(c(n_timepoints = 0, sd_subject = NA_real_, range_subject = NA_real_, mssd = NA_real_))
    }
    sdv <- if (length(x) >= 2) sd(x) else NA_real_
    rgv <- if (length(x) >= 1) max(x) - min(x) else NA_real_
    msv <- if (length(x) >= 2) mean(diff(x)^2) else NA_real_
    c(n_timepoints = length(x), sd_subject = sdv, range_subject = rgv, mssd = msv)
  }
  
  # One row per participant × measure with SD, Range, MSSD
  subj_tbl <- long_df %>%
    dplyr::filter(.data$Measure %in% measures) %>%
    dplyr::arrange(.data[[id_col]], .data$Measure, .data[[time_col]]) %>%
    dplyr::group_by(.data[[id_col]], .data[[group_col]], .data$Measure) %>%
    dplyr::summarise(metrics = list(subject_metrics(.data[[value_col]])),
                     .groups = "drop") %>%
    tidyr::unnest_wider(metrics)
  
  # Wilcoxon tests (cases vs controls) for SD, Range, MSSD per measure
  metrics_to_test <- c("sd_subject", "range_subject", "mssd")
  wilcox_rows <- list()
  
  for (m in unique(subj_tbl$Measure)) {
    dat_m <- subj_tbl %>% dplyr::filter(Measure == m)
    if (dplyr::n_distinct(dat_m[[group_col]]) < 2) next
    
    for (met in metrics_to_test) {
      d2 <- dat_m %>% dplyr::filter(!is.na(.data[[met]]))
      
      if (dplyr::n_distinct(d2[[group_col]]) < 2) {
        p <- NA_real_; W <- NA_real_
      } else {
        test <- wilcox.test(
          stats::as.formula(paste(met, "~", group_col)),
          data = d2,
          exact = exact_two_sided
        )
        p <- unname(test$p.value)
        W <- unname(test$statistic)
      }
      
      wilcox_rows[[length(wilcox_rows) + 1]] <- tibble::tibble(
        Measure = m,
        Metric = met,
        W = W,
        p_value = p,
        n_subjects_total = nrow(d2),
        n_groups = dplyr::n_distinct(d2[[group_col]])
      )
    }
  }
  
  wilcox_tbl <- dplyr::bind_rows(wilcox_rows) %>%
    dplyr::mutate(Metric = dplyr::recode(Metric,
                                         sd_subject = "SD",
                                         range_subject = "Range",
                                         mssd = "MSSD"))
  
  list(subject_metrics = subj_tbl, wilcox_tests = wilcox_tbl)
}


#Step 2b: run
TMT_measures <- c("tmt_trailA_reactiontime", "tmt_trailAminustrailB_exefcn_final", "tmt_trailB_reactiontime")

res_var <- calc_subject_variability(
  long_df   = sRBD_hyp_TMT_long_all,
  measures  = TMT_measures,
  value_col = "Value",
  exact_two_sided = FALSE
)

subject_variance_all <- res_var$subject_metrics
stat_results_all     <- res_var$wilcox_tests

#save
write.csv(subject_variance_all,
          "TMT/sRBD hyp vs controls/subject_variance_all_measures_sRBD_hyp.csv",
          row.names = FALSE)

write.csv(stat_results_all,
          "TMT/sRBD hyp vs controls/variability_stats_wilcoxon_all_measures_sRBD_hyp.csv",
          row.names = FALSE)

#create boxplots
make_variability_boxplot <- function(subject_tbl, measure_name, plot_title) {
  subject_tbl %>%
    filter(Measure == measure_name) %>%
    pivot_longer(
      cols = c(sd_subject, range_subject, mssd),
      names_to = "Metric",
      values_to = "Value_plot"
    ) %>%
    mutate(
      Metric = recode(Metric,
                      sd_subject = "SD",
                      range_subject = "Range",
                      mssd = "MSSD")
    ) %>%
    ggplot(aes(x = Code, y = Value_plot, fill = Code)) +
    geom_boxplot(alpha = 0.5) +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(title = plot_title, x = "Group", y = "Value") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
}


p_ART <- make_variability_boxplot(subject_variance_all,
                                  "tmt_trailA_reactiontime",
                                  "Within-Subject Variability of TMT-A RT")

p_EF <- make_variability_boxplot(subject_variance_all,
                                 "tmt_trailAminustrailB_exefcn_final",
                                 "Within-Subject Variability of Executive Function")

p_BRT <- make_variability_boxplot(subject_variance_all,
                                  "tmt_trailB_reactiontime",
                                  "Within-Subject Variability of TMT-B RT")

print(p_ART); print(p_EF); print(p_BRT)
subject_variance_all %>% count(Measure, n_timepoints)

wilcox_table <- stat_results_all %>%
  select(Measure, Metric, p_value) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = p_value)

wilcox_table


#Step 3:	Model time-varying variance (heteroscedastic mixed model)
#Filter Trail A RT scores
taRT_data <- sRBD_hyp_TMT_long_all %>% 
  filter(Measure == "tmt_trailA_reactiontime")
#make relative visit a factor 
taRT_data <- taRT_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))


#Homoscedastic LME (residual variance constant)
lme_hom <- nlme::lme(
  Value ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = taRT_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident <- update(
  lme_hom,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom, lme_varident)

summary(lme_hom)


# Filter Trail B RT scores
tbRT_data <- sRBD_hyp_TMT_long_all %>% 
  filter(Measure == "tmt_trailB_reactiontime")
#make relative visit a factor 
tbRT_data <- tbRT_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))

#Homoscedastic LME (residual variance constant)
lme_hom_bRT <- nlme::lme(
  Value ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = tbRT_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident_bRT <- update(
  lme_hom_bRT,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom_bRT, lme_varident_bRT)

summary(lme_varident_bRT)


# Filter EF scores
EF_data <- sRBD_hyp_TMT_long_all %>% 
  filter(Measure == "tmt_trailAminustrailB_exefcn_final")
#make relative visit a factor 
EF_data <- EF_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))

#Homoscedastic LME (residual variance constant)
lme_hom_EF <- nlme::lme(
  Value ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = EF_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident_EF <- update(
  lme_hom_EF,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom_EF, lme_varident_EF)

summary(lme_varident_EF)


###### demographics ###### 
sRBD_hyp_TMT_long_all  %>%
  distinct(PID2, Code) %>%
  count(Code)

#keep one row per participant for demographics
demo_df <- sRBD_hyp_TMT_long_all %>%
  arrange(PID2, RelativeVisit) %>% 
  distinct(PID2, .keep_all = TRUE)

t_test_age_result <- t.test(age_session_baseline ~ Code, data = demo_df)
print(t_test_age_result)

std_dev_ages <- demo_df %>% group_by(Code) %>% 
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE))
print(std_dev_ages)

sex_percentage <- demo_df %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

chi_test_result_sex <- demo_df %>%
  with(table(gender, Code)) %>%
  chisq.test()
print(chi_test_result_sex)

summary_means <- sRBD_hyp_TMT_long_all %>%
  group_by(Measure, RelativeVisit) %>%
  summarise(
    n_control = sum(Code == "Control" & !is.na(Value)),
    mean_control = mean(Value[Code == "Control"], na.rm = TRUE),
    sd_control = sd(Value[Code == "Control"], na.rm = TRUE),
    n_sRBD_hyp = sum(Code == "sRBD_hyposmia" & !is.na(Value)),
    mean_sRBD_hyp = mean(Value[Code == "sRBD_hyposmia"], na.rm = TRUE),
    sd_sRBD_hyp = sd(Value[Code == "sRBD_hyposmia"], na.rm = TRUE),
    p_value = wilcox.test(Value ~ Code, data = pick(everything()), exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    control_mean_sd = sprintf("%.2f ± %.2f", mean_control, sd_control),
    sRBD_hyp_mean_sd = sprintf("%.2f ± %.2f", mean_sRBD_hyp, sd_sRBD_hyp)
  ) %>%
  select(Measure, RelativeVisit, n_control, control_mean_sd,
         n_sRBD_hyp, sRBD_hyp_mean_sd, p_value)

print(summary_means, n = nrow(summary_means))
print(as.data.frame(summary_means), digits = 10)

######### DFT sRBD vs controls baseline analysis ##########
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(moments)
library(logistf)
library(broom) # to tidy
library(ResourceSelection) # for Hosmer-Lemeshow test 

sRBD_control_DFT_data <- read_csv("DFT/sRBD_control_DFT_data.csv")

table(sRBD_control_DFT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
sRBD_control_DFT_data <- sRBD_control_DFT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))

## remove those with essential tremor and under 50
sRBD_control_DFT_data <- sRBD_control_DFT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)


sRBD_control_DFT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)
# remove high risk controls
sRBD_control_DFT_data <- sRBD_control_DFT_data %>%
  filter(!(Code == "Control" & PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline == "Yes"))
# double check
sRBD_control_DFT_data %>%
  filter(Code == "Control") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)


#get rid of IDs with missing data
sRBD_control_baseline_data <- sRBD_control_DFT_data %>%
  filter(!is.na(onefingertaptest_worst_ks_baseline) & !is.na(onefingertaptest_worst_is_baseline) & !is.na(onefingertaptest_worst_at_baseline))


table(sRBD_control_baseline_data$Code)

## reverse direction of KS
sRBD_control_baseline_data$onefingertaptest_worst_ks_rev <- -sRBD_control_baseline_data$onefingertaptest_worst_ks_baseline
#check 
head(sRBD_control_baseline_data[, c("onefingertaptest_worst_ks_baseline", "onefingertaptest_worst_ks_rev")])


### make subgroup for those with PD
sRBD_and_PD_bl_data <- sRBD_control_baseline_data %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
sRBD_control_baseline_data <- sRBD_control_baseline_data %>%
  anti_join(sRBD_and_PD_bl_data, by = c("PID2", "TokenFinal"))
table(sRBD_control_baseline_data$hx_neuro_parkinsons, sRBD_control_baseline_data$Code)

write.csv(sRBD_control_baseline_data, "DFT/Pilot study dft/Baseline/sRBD vs controls/sRBD_control_baseline_data.csv", row.names = FALSE)




#### baseline analysis sRBD vs controls ####
summary(sRBD_control_baseline_data)
#get rid of IDs with missing data
sRBD_control_baseline_data <- sRBD_control_baseline_data %>%
  filter(!is.na(onefingertaptest_worst_ks_rev) & !is.na(onefingertaptest_worst_is_baseline) & !is.na(onefingertaptest_worst_at_baseline))
#get rid of '0'
sRBD_control_baseline_data <- sRBD_control_baseline_data %>%
  filter(
    abs(onefingertaptest_worst_is_baseline) > 1e-6,
    abs(onefingertaptest_worst_at_baseline) > 1e-6,
    abs(onefingertaptest_worst_ks_rev) > 1e-6
  )

table(sRBD_control_baseline_data$Code)


# log transform IS
sRBD_control_baseline_data$onefingertaptest_worst_is_baseline <- log(sRBD_control_baseline_data$onefingertaptest_worst_is_baseline)
#remove NaN rows
sRBD_control_baseline_data2 <- sRBD_control_baseline_data %>%
  filter(!is.nan(onefingertaptest_worst_is_baseline))

#have outliers in control group but don't want to 'stack the deck' to get sig differences
#instead of excluding set the value at the cohort max/min
# Winsorization function to replace outliers with the max/min within bounds
winsorize_with_max_min <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lower <- qnt[1] - 1.5 * iqr
  upper <- qnt[2] + 1.5 * iqr
  
  # Find the maximum and minimum values within the bounds
  capped_min <- min(x[x >= lower], na.rm = TRUE)
  capped_max <- max(x[x <= upper], na.rm = TRUE)
  
  # Replace outliers with the max/min within the bounds
  x[x < lower] <- capped_min
  x[x > upper] <- capped_max
  
  return(x)
}

# Apply to the selected columns to controls
sRBD_control_baseline_data2 <- sRBD_control_baseline_data2 %>%
  group_by(Code) %>%
  mutate(across(
    c(onefingertaptest_worst_at_baseline,
      onefingertaptest_worst_is_baseline,
      onefingertaptest_worst_ks_rev),
    ~ {
      if (first(Code) == "Control") winsorize_with_max_min(.) else .
    }
  )) %>%
  ungroup()


# check categories of data
str(sRBD_control_baseline_data2)
table(sRBD_control_baseline_data2$Code)



#### check distribution of baseline data - density plots ####
ggplot(sRBD_control_baseline_data2, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45, 87)


DFT_baseline_male <- sRBD_control_baseline_data2[sRBD_control_baseline_data2$gender == "Male", ]

ggplot(DFT_baseline_male, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,88)


DFT_baseline_female <- sRBD_control_baseline_data2[sRBD_control_baseline_data2$gender == "Female", ]

ggplot(DFT_baseline_female, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,88)

ggplot(sRBD_control_baseline_data2, aes(x=onefingertaptest_worst_is_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(4, 13)

ggplot(sRBD_control_baseline_data2, aes(x=onefingertaptest_worst_at_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(0, 300)

ggplot(sRBD_control_baseline_data2, aes(x=onefingertaptest_worst_ks_rev, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(-200, 25)




#### Demographics ####
table(sRBD_control_baseline_data2$Code)
t_test_age_result <- t.test(age_session_baseline ~ Code, data = sRBD_control_baseline_data2)
print(t_test_age_result)

std_dev_ages <- sRBD_control_baseline_data2 %>%
  group_by(Code) %>%
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE)) %>%
  ungroup()
print(std_dev_ages)

sex_percentage <- sRBD_control_baseline_data2 %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

contingency_table <- table(sRBD_control_baseline_data2$gender, sRBD_control_baseline_data2$Code)
chi_test_result_sex <- chisq.test(contingency_table)
print(chi_test_result_sex)

sRBD_control_baseline_data2 %>%
  group_by(Code) %>%
  summarise(median_score = median(onefingertaptest_worst_is_baseline, na.rm = TRUE), .groups = "drop")
wilcox.test(onefingertaptest_worst_is_baseline ~ Code, data = sRBD_control_baseline_data2)

t.test(onefingertaptest_worst_at_baseline ~ Code, data = sRBD_control_baseline_data2)

t.test(onefingertaptest_worst_ks_rev ~ Code, data = sRBD_control_baseline_data2)

#PD
table(sRBD_control_baseline_data2$Code, sRBD_control_baseline_data2$hx_neuro_parkinsons)



#### Regression model ####
# convert code column to 0 (control) and 1 (iRBD) for glm
sRBD_control_baseline_data2 <- sRBD_control_baseline_data2 %>%
  mutate(Code = ifelse(Code == "sRBD", 1, 0))
sRBD_control_baseline_data2$Code <- as.factor(sRBD_control_baseline_data2$Code)
table(sRBD_control_baseline_data2$Code)


# Incoordination score
ggplot(sRBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_is_baseline)) + 
  geom_boxplot()
# standardise and use z score
sRBD_control_baseline_data2$onefingertaptest_worst_is_baseline_z <- scale(sRBD_control_baseline_data2$onefingertaptest_worst_is_baseline)
model_IS <- glm( Code ~ onefingertaptest_worst_is_baseline_z  + gender + age_session_baseline, 
                 data = sRBD_control_baseline_data2,family = 'binomial')
summary(model_IS)
exp(cbind(coef(model_IS), confint(model_IS,level = 0.95)))



# Kinesia score
ggplot(sRBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_ks_rev)) + 
  geom_boxplot()
# standardise and use z score
sRBD_control_baseline_data2$onefingertaptest_worst_ks_rev_z <- scale(sRBD_control_baseline_data2$onefingertaptest_worst_ks_rev)
model_KS <- glm( Code ~ onefingertaptest_worst_ks_rev_z  + gender + age_session_baseline, 
                 data = sRBD_control_baseline_data2,family = 'binomial')
summary(model_KS)
exp(cbind(coef(model_KS), confint(model_KS,level = 0.95)))



#Akinesia time
ggplot(sRBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_at_baseline)) + 
  geom_boxplot()
# standardise and use z score
sRBD_control_baseline_data2$onefingertaptest_worst_at_baseline_z <- scale(sRBD_control_baseline_data2$onefingertaptest_worst_at_baseline)
model_AT <- glm( Code ~ onefingertaptest_worst_at_baseline_z  + gender + age_session_baseline, 
                 data = sRBD_control_baseline_data2,family = 'binomial')
summary(model_AT)
exp(cbind(coef(model_AT), confint(model_AT,level = 0.95)))


### check age linearity
#IS
ggplot(sRBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_IS))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, IS", x = "Age", y = "Predicted Probability of Code = 1")

sRBD_control_baseline_data2$age_baseline_sq <- sRBD_control_baseline_data2$age_session_baseline^2

model_age_quad <- glm(Code ~ onefingertaptest_worst_is_baseline_z + gender + age_session_baseline + age_baseline_sq, 
                      data = sRBD_control_baseline_data2, family = "binomial")
summary(model_age_quad)
exp(cbind(coef(model_age_quad), confint(model_age_quad,level = 0.95)))
AIC(model_IS, model_age_quad)
final_model_IS <- model_age_quad


#KS
ggplot(sRBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_KS))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, KS", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_KS <- glm(Code ~ onefingertaptest_worst_ks_rev_z + gender + age_session_baseline + age_baseline_sq, 
                         data = sRBD_control_baseline_data2, family = "binomial")
summary(model_age_quad_KS)
exp(cbind(coef(model_age_quad_KS), confint(model_age_quad_KS,level = 0.95)))
AIC(model_KS, model_age_quad_KS)
final_model_KS <- model_age_quad_KS



#AT
ggplot(sRBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_AT))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, AT", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_AT <- glm(Code ~ onefingertaptest_worst_at_baseline_z + gender + age_session_baseline + age_baseline_sq, 
                         data = sRBD_control_baseline_data2, family = "binomial")
summary(model_age_quad_AT)
exp(cbind(coef(model_age_quad_AT), confint(model_age_quad_AT,level = 0.95)))
AIC(model_AT, model_age_quad_AT)
final_model_AT <- model_age_quad_AT


## QQ plots ##
# List of models and labels
models <- list(
  IS = final_model_IS,
  KS = final_model_KS,
  AT = final_model_AT
)


df_all_residuals <- bind_rows(
  lapply(names(models), function(name) {
    data.frame(
      residuals = residuals(models[[name]], type = "deviance"),
      Outcome = name,
      Code = sRBD_control_baseline_data2$Code
    )
  })
)

ggplot(df_all_residuals, aes(sample = residuals, color = factor(Code))) +
  stat_qq() +
  stat_qq_line(color = "black") +
  facet_wrap(~ Outcome, scales = "free") +
  labs(
    title = "QQ Plots of Deviance Residuals by Outcome",
    color = "Group",
    x = "Theoretical Quantiles",
    y = "Deviance Residuals"    
  ) +
  theme_minimal()



### tidy up final models
tidy(final_model_IS, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)

tidy(final_model_KS, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)

tidy(final_model_AT, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)




#### Hosmer-Lemeshow test to check if the model fits ####
#IS
sRBD_control_baseline_data2$Code_numeric <- as.numeric(as.character(sRBD_control_baseline_data2$Code))

hoslem.test(x = sRBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_IS),
            g = 10)

#KS
hoslem.test(x = sRBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_KS),
            g = 10)


#AT
hoslem.test(x = sRBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_AT),
            g = 10)





##### Create Forest plots ##### 
#IS
# Extract ORs and 95% CIs
OR_data <- as.data.frame(exp(cbind(OR = coef(final_model_IS), confint(final_model_IS))))
OR_data$Variable <- rownames(OR_data)

# Replace raw variable names with  labels
OR_data$Variable <- recode(OR_data$Variable,
                           "(Intercept)" = "Intercept",
                           "onefingertaptest_worst_is_baseline_z" = "IS (z-score)",
                           "genderMale" = "Male (vs Female)",
                           "age_session_baseline" = "Age",
                           "age_baseline_sq" = "Age²"
)



# Exclude Intercept for clarity
ggplot(OR_data[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, linewidth = 1.2, # <- increase size for boldness
                color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +  # <- optional: make dashed line bolder
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)



#KS
OR_data_KS <- as.data.frame(exp(cbind(OR = coef(final_model_KS), confint(final_model_KS))))
OR_data_KS$Variable <- rownames(OR_data_KS)

# Replace raw variable names with  labels
OR_data_KS$Variable <- recode(OR_data_KS$Variable,
                              "(Intercept)" = "Intercept",
                              "onefingertaptest_worst_ks_rev_z" = "KS (z-score)",
                              "genderMale" = "Male (vs Female)",
                              "age_session_baseline" = "Age",
                              "age_baseline_sq" = "Age²"
)


ggplot(OR_data_KS[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, size = 1.2,
                color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +  
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)



#AT
OR_data_AT <- as.data.frame(exp(cbind(OR = coef(final_model_AT), confint(final_model_AT))))
OR_data_AT$Variable <- rownames(OR_data_AT)

# Replace raw variable names with  labels
OR_data_AT$Variable <- recode(OR_data_AT$Variable,
                              "(Intercept)" = "Intercept",
                              "onefingertaptest_worst_at_baseline_z" = "AT (z-score)",
                              "genderMale" = "Male (vs Female)",
                              "age_session_baseline" = "Age",
                              "age_baseline_sq" = "Age²"
)

ggplot(OR_data_AT[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, size = 1.2,
                color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +  
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)




# collate into one forest plot
ks_data <- tidy(final_model_KS, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
  filter(term == "onefingertaptest_worst_ks_rev_z") %>%
  mutate(Variable = "Kinesia Score (KS)")

is_data <- tidy(final_model_IS, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
  filter(term == "onefingertaptest_worst_is_baseline_z") %>%
  mutate(Variable = "Incoordination Score (IS)")

at_data <- tidy(final_model_AT, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
  filter(term == "onefingertaptest_worst_at_baseline_z") %>%
  mutate(Variable = "Akinesia Time (AT)")

all_motor <- bind_rows(ks_data, is_data, at_data) %>%
  mutate(p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))



ggplot(all_motor, aes(x = estimate, y = Variable)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "#2c7fb8") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 1) +
  geom_text(aes(label = paste("p =", p_label), x = conf.high * 1.15), 
            hjust = 0, size = 4.2, color = "black") +
  scale_x_log10(limits = c(0.50, 4)) +  # ← sets fixed axis range
  labs(title = "Baseline Associations of DFT Scores with normosmic pRBD vs. Controls",
       subtitle = "Odds Ratios from Logistic Regression Adjusted for Age and Gender",
       x = "Odds Ratio (log scale)", y = "") +
  theme_minimal(base_size = 14)




#extract data for overall forest plot
extract_row <- function(model, term, variable, test, comparison,
                        conf_level = 0.95) {
  tidy(model, conf.int = TRUE, conf.level = conf_level, exponentiate = TRUE) %>%
    filter(term == !!term) %>%
    mutate(
      Variable   = variable,
      Test       = test,
      Comparison = comparison
    ) %>%
    select(Test, Comparison, Variable, estimate, conf.low, conf.high, p.value)
}


DFT_df_sRBD <- bind_rows(
  extract_row(final_model_KS,
              "onefingertaptest_worst_ks_rev_z",
              "Kinesia Score",
              test = "DFT", comparison = "Normosmic pRBD vs control"),
  extract_row(final_model_IS,
              "onefingertaptest_worst_is_baseline_z",
              "Incoordination Score",
              test = "DFT", comparison = "Normosmic pRBD vs control"),
  extract_row(final_model_AT,
              "onefingertaptest_worst_at_baseline_z",
              "Akinesia Time",
              test = "DFT", comparison = "Normosmic pRBD vs control")
)

saveRDS(DFT_df_sRBD, "DFT/Pilot study dft/Baseline/sRBD vs controls/forest_rows_DFT_sRBD.rds")












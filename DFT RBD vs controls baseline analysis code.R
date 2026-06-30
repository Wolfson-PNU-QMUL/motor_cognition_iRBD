######### DFT RBD vs controls baseline analysis ##########
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(moments)
library(broom) 
library(ResourceSelection)  

table(RBD_control_DFT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))


## remove those with essential tremor and under 50
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)


RBD_control_DFT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)
# remove high risk controls
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  filter(!(Code == "Control" & PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline == "Yes"))
# double check
RBD_control_DFT_data %>%
  filter(Code == "Control") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)

#get rid of IDs with missing data
RBD_control_baseline_data <- RBD_control_DFT_data %>%
  filter(!is.na(onefingertaptest_worst_ks_baseline) & !is.na(onefingertaptest_worst_is_baseline) & !is.na(onefingertaptest_worst_at_baseline))


table(RBD_control_baseline_data$Code)

## reverse direction of KS
RBD_control_baseline_data$onefingertaptest_worst_ks_rev <- -RBD_control_baseline_data$onefingertaptest_worst_ks_baseline
#check 
head(RBD_control_baseline_data[, c("onefingertaptest_worst_ks_baseline", "onefingertaptest_worst_ks_rev")])

### make subgroup for those with PD
RBD_and_PD_bl_data <- RBD_control_baseline_data %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
RBD_control_baseline_data <- RBD_control_baseline_data %>%
  anti_join(RBD_and_PD_bl_data, by = c("PID2", "TokenFinal"))



#### baseline analysis RBD vs controls ####
summary(RBD_control_baseline_data)
#get rid of IDs with missing data
RBD_control_baseline_data <- RBD_control_baseline_data %>%
  filter(!is.na(onefingertaptest_worst_ks_rev) & !is.na(onefingertaptest_worst_is_baseline) & !is.na(onefingertaptest_worst_at_baseline))
#get rid of '0'
RBD_control_baseline_data <- RBD_control_baseline_data %>%
  filter(
    abs(onefingertaptest_worst_is_baseline) > 1e-6,
    abs(onefingertaptest_worst_at_baseline) > 1e-6,
    abs(onefingertaptest_worst_ks_rev) > 1e-6
  )

table(RBD_control_baseline_data$Code)

# log transform IS
RBD_control_baseline_data$onefingertaptest_worst_is_baseline <- log(RBD_control_baseline_data$onefingertaptest_worst_is_baseline)
#remove NaN rows
RBD_control_baseline_data2 <- RBD_control_baseline_data %>%
  filter(!is.nan(onefingertaptest_worst_is_baseline))


#have outliers in control group
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
RBD_control_baseline_data2 <- RBD_control_baseline_data2 %>%
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
str(RBD_control_baseline_data2)
table(RBD_control_baseline_data2$Code)

#### check distribution of baseline data - density plots ####
ggplot(RBD_control_baseline_data2, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,90)

DFT_baseline_male <- RBD_control_baseline_data2[RBD_control_baseline_data2$gender == "Male", ]

ggplot(DFT_baseline_male, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,93)

DFT_baseline_female <- RBD_control_baseline_data2[RBD_control_baseline_data2$gender == "Female", ]

ggplot(DFT_baseline_female, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,87)


ggplot(RBD_control_baseline_data2, aes(x=onefingertaptest_worst_is_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(4, 12)
IS <- RBD_control_baseline_data2$onefingertaptest_worst_is_baseline
IS <- rnorm(5802)
skew_value <- skewness(IS)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}

ggplot(RBD_control_baseline_data2, aes(x=onefingertaptest_worst_at_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(0, 275)

ggplot(RBD_control_baseline_data2, aes(x=onefingertaptest_worst_ks_rev, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(-150, 30)



#### Demographics ####
t_test_age_result <- t.test(age_session_baseline ~ Code, data = RBD_control_baseline_data2)
print(t_test_age_result)

std_dev_ages <- RBD_control_baseline_data2 %>%
  group_by(Code) %>%
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE)) %>%
  ungroup()
print(std_dev_ages)

sex_percentage <- RBD_control_baseline_data2 %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

contingency_table <- table(RBD_control_baseline_data2$gender, RBD_control_baseline_data2$Code)
chi_test_result_sex <- chisq.test(contingency_table)
print(chi_test_result_sex)

RBD_control_baseline_data2 %>%
  group_by(Code) %>%
  summarise(median_score = median(onefingertaptest_worst_is_baseline, na.rm = TRUE), .groups = "drop")
wilcox.test(onefingertaptest_worst_is_baseline ~ Code, data = RBD_control_baseline_data2)

t.test(onefingertaptest_worst_at_baseline ~ Code, data = RBD_control_baseline_data2)

t.test(onefingertaptest_worst_ks_rev ~ Code, data = RBD_control_baseline_data2)

#PD
table(RBD_control_baseline_data2$Code, RBD_control_baseline_data2$hx_neuro_parkinsons)
PD_table <- table(RBD_control_baseline_data2$Code, RBD_control_baseline_data2$hx_neuro_parkinsons)
# Fisher’s Exact Test
fisher.test(PD_table)

#### Regression model ####
# convert code column to 0 (control) and 1 (iRBD) for glm
RBD_control_baseline_data2 <- RBD_control_baseline_data2 %>%
  mutate(Code = ifelse(Code == "iRBD", 1, 0))
RBD_control_baseline_data2$Code <- as.factor(RBD_control_baseline_data2$Code)
table(RBD_control_baseline_data2$Code)


# Incoordination score
ggplot(RBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_is_baseline)) + 
  geom_boxplot()
# standardise and use z score
RBD_control_baseline_data2$onefingertaptest_worst_is_baseline_z <- scale(RBD_control_baseline_data2$onefingertaptest_worst_is_baseline)
model_IS <- glm( Code ~ onefingertaptest_worst_is_baseline_z  + gender + age_session_baseline, 
                 data = RBD_control_baseline_data2,family = 'binomial')
summary(model_IS)
exp(cbind(coef(model_IS), confint(model_IS,level = 0.95)))


# Kinesia score
ggplot(RBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_ks_rev)) + 
  geom_boxplot()
# standardise and use z score
RBD_control_baseline_data2$onefingertaptest_worst_ks_rev_z <- scale(RBD_control_baseline_data2$onefingertaptest_worst_ks_rev)
model_KS <- glm( Code ~ onefingertaptest_worst_ks_rev_z  + gender + age_session_baseline, 
                 data = RBD_control_baseline_data2,family = 'binomial')
summary(model_KS)
exp(cbind(coef(model_KS), confint(model_KS,level = 0.95)))


#Akinesia time
ggplot(RBD_control_baseline_data2, aes(x=Code, y=onefingertaptest_worst_at_baseline)) + 
  geom_boxplot()
# standardise and use z score
RBD_control_baseline_data2$onefingertaptest_worst_at_baseline_z <- scale(RBD_control_baseline_data2$onefingertaptest_worst_at_baseline)
model_AT <- glm( Code ~ onefingertaptest_worst_at_baseline_z  + gender + age_session_baseline, 
                 data = RBD_control_baseline_data2,family = 'binomial')
summary(model_AT)
exp(cbind(coef(model_AT), confint(model_AT,level = 0.95)))



### check age linearity
#IS
ggplot(RBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_IS))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, IS", x = "Age", y = "Predicted Probability of Code = 1")

RBD_control_baseline_data2$age_baseline_sq <- RBD_control_baseline_data2$age_session_baseline^2

model_age_quad <- glm(Code ~ onefingertaptest_worst_is_baseline_z + gender + age_session_baseline + age_baseline_sq, 
                      data = RBD_control_baseline_data2, family = "binomial")
summary(model_age_quad)
exp(cbind(coef(model_age_quad), confint(model_age_quad,level = 0.95)))
AIC(model_IS, model_age_quad)
final_model_IS <- model_age_quad


#KS
ggplot(RBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_KS))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, KS", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_KS <- glm(Code ~ onefingertaptest_worst_ks_rev_z + gender + age_session_baseline + age_baseline_sq, 
                         data = RBD_control_baseline_data2, family = "binomial")
summary(model_age_quad_KS)
exp(cbind(coef(model_age_quad_KS), confint(model_age_quad_KS,level = 0.95)))
AIC(model_KS, model_age_quad_KS)
final_model_KS <- model_age_quad_KS



#AT
ggplot(RBD_control_baseline_data2, aes(x = age_session_baseline, y = fitted(model_AT))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, AT", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_AT <- glm(Code ~ onefingertaptest_worst_at_baseline_z + gender + age_session_baseline + age_baseline_sq, 
                         data = RBD_control_baseline_data2, family = "binomial")
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
      Code = RBD_control_baseline_data2$Code
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
RBD_control_baseline_data2$Code_numeric <- as.numeric(as.character(RBD_control_baseline_data2$Code))

hoslem.test(x = RBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_IS),
            g = 10)

#KS
hoslem.test(x = RBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_KS),
            g = 10)

#AT
hoslem.test(x = RBD_control_baseline_data2$Code_numeric,
            y = fitted(final_model_AT),
            g = 10)

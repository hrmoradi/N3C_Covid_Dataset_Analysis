library(MatchIt)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

@transform_pandas(
    Output(rid="ri.vector.main.execute.cf865dec-7d73-4ac1-b16c-d2aa15d58fcb"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
Propensity_Matching_GFR_Impute <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_imputed = ifelse(is.na(GFR), 1, 0))
momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    Num_GFR = (n() - sum(GFR_imputed)),
    Num_GFR_imputed = sum(GFR_imputed)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR','GFR_imputed')

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR +GFR_imputed,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR + GFR_imputed,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' , 'GFR_imputed','distance')
diff_means <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
   fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
   fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 4, widths = c(1, 0.85)
)
# return() ################## 1

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR+ GFR_imputed, family = binomial, data = dta_m)
print(summary(glm_treat2))

# plotting Residuals
print("*** plotting Residuals")
plot(predict(glm_treat2),residuals(glm_treat2))

# Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.677e25ca-82a5-43b1-b7f6-5d1d3ed7a0f1"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
Propensity_Matching_overal <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
#momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients))#,
    # sum_GFR_exist = ( sum(GFR_exist)),
    # sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ) # ,'GFR'

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  , # + GFR 
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   , # + GFR
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

tab_after <- map_df(diff_means_after, tidy)

tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',  'distance' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
   #fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
#    fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 3, widths = c(1, 0.85)
)
# return() ##################
return(tab)

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.e4ee75aa-d729-4bca-a395-5a8c9b5ed999"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
Propensity_Matching_overal_GLM <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
#momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients))#,
    # sum_GFR_exist = ( sum(GFR_exist)),
    # sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ) # ,'GFR'

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  , # + GFR 
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   , # + GFR
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))
return(dta_m)
# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# tab_after <- map_df(diff_means_after, tidy)

# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',  'distance' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# Checking balance
library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    #fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# return() ##################
# return(tab)

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7fc019ac-b7c6-4e1a-84f8-4589f8196eec"),
    results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.dc40e44c-80a4-4e85-86bf-07b4a9ed206a")
)
after_match_Cohort_negative_with_GFR_treatment_Effect_ttest <- function(results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('mort_hospice') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('mort_hospice'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # # plotting Residuals
# # print("*** plotting Residuals")
# # plot(predict(glm_treat2),residuals(glm_treat2))

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.122f0b9e-6767-4240-af4f-cdcb3ec10bc4"),
    results_Propensity_Matching_cohort_without_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.77db6ef3-28d6-4ce8-8ef4-536367d5028d")
)
after_match_Cohort_negative_without_GFR_treatment_Effect_ttest <- function(results_Propensity_Matching_cohort_without_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_LogisticModel

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('mort_hospice') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('mort_hospice'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # # plotting Residuals
# # print("*** plotting Residuals")
# # plot(predict(glm_treat2),residuals(glm_treat2))

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.0c666b90-557b-4e27-9617-9fcf8fcfa3de"),
    results_Propensity_Matching_cohort_with_GFR_GLM=Input(rid="ri.foundry.main.dataset.00394da4-c9bc-4a9d-9776-db29ca7c0464")
)
after_match_Cohort_with_GFR_treatment_Effect <- function(results_Propensity_Matching_cohort_with_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_with_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('mort_hospice') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('mort_hospice'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # # plotting Residuals
# # print("*** plotting Residuals")
# # plot(predict(glm_treat2),residuals(glm_treat2))


# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)


}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.8a215568-2de6-4a1d-8f3d-292f3fb716b4"),
    results_Propensity_Matching_cohort_without_GFR_GLM=Input(rid="ri.foundry.main.dataset.7bbd5171-673f-4bf0-8522-a9f3419123f8")
)
after_match_Cohort_without_GFR_treatment_Effect <- function(results_Propensity_Matching_cohort_without_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('mort_hospice') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('mort_hospice'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # # plotting Residuals
# # print("*** plotting Residuals")
# # plot(predict(glm_treat2),residuals(glm_treat2))

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.51fd7c99-0d95-4530-931e-7b89b1c91ad8"),
    Propensity_Matching_overal_GLM=Input(rid="ri.foundry.main.dataset.e4ee75aa-d729-4bca-a395-5a8c9b5ed999")
)
after_match_overal_GLM_ACE_ARB_only <- function(Propensity_Matching_overal_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- Propensity_Matching_overal_GLM

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.1b3c7d2e-7cfc-4b37-8ca7-59762f8c6ad7"),
    Propensity_Matching_overal_GLM=Input(rid="ri.foundry.main.dataset.e4ee75aa-d729-4bca-a395-5a8c9b5ed999")
)
after_match_overal_GLM_ALL <- function(Propensity_Matching_overal_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- Propensity_Matching_overal_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)
# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.35bf8c99-0430-4094-b3b1-2cf24bfc1ede"),
    Propensity_Matching_overal_GLM=Input(rid="ri.foundry.main.dataset.e4ee75aa-d729-4bca-a395-5a8c9b5ed999")
)
after_match_overal_GLM_sudo_r <- function(Propensity_Matching_overal_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- Propensity_Matching_overal_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm - only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

Vars <- c("insample glm - only ACE/ARB", "insample glm all variables ")
Pseudo_R_Squared <- c(1 - glm_treat1$deviance / glm_treat1$null.deviance, 1 - glm_treat2$deviance / glm_treat2$null.deviance)

df <- data.frame(Vars, Pseudo_R_Squared)
return(df)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f23e2490-c284-40bc-961c-3dffbe1f1bca"),
    Propensity_Matching_overal_GLM=Input(rid="ri.foundry.main.dataset.e4ee75aa-d729-4bca-a395-5a8c9b5ed999")
)
after_match_overal_treatment_Effect <- function(Propensity_Matching_overal_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- Propensity_Matching_overal_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('mort_hospice') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('mort_hospice'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

# # # plotting Residuals
# # print("*** plotting Residuals")
# # plot(predict(glm_treat2),residuals(glm_treat2))

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c52d98f5-46f5-4f30-ae47-7a4aa5495ee2"),
    results_Propensity_Matching_cohort_with_GFR_GLM=Input(rid="ri.foundry.main.dataset.00394da4-c9bc-4a9d-9776-db29ca7c0464")
)
after_match_with_GFR_GLM_ACE_ARB_only <- function(results_Propensity_Matching_cohort_with_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_with_GFR_GLM

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.9755ed0f-2123-45d6-ae4b-974de4bdd546"),
    results_Propensity_Matching_cohort_with_GFR_GLM=Input(rid="ri.foundry.main.dataset.00394da4-c9bc-4a9d-9776-db29ca7c0464")
)
after_match_with_GFR_GLM_ALL <- function(results_Propensity_Matching_cohort_with_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_with_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR   , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)
# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.81d59d77-7a49-49fc-9f9e-ca5650180b0b"),
    results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.dc40e44c-80a4-4e85-86bf-07b4a9ed206a")
)
after_match_with_GFR_LogisticModel_ACE_ARB_only <- function(results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.156a38b2-fc08-4603-8cfd-8d5256f95c12"),
    results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.dc40e44c-80a4-4e85-86bf-07b4a9ed206a")
)
after_match_with_GFR_LogisticModel_ALL <- function(results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# tab <- tidy(glm_treat1)
# tab <- tab %>% 
#   rename(
#     # mean_difference = estimate,
#     z_statistic = statistic
#     # df_parameter = parameter,
#     # mean_not_ACE_ARB = estimate1,
#     # Pr_z = p_value
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR   , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)
# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.cca38d94-9419-4cbc-b353-8e28d858546b"),
    results_Propensity_Matching_cohort_with_GFR_GLM=Input(rid="ri.foundry.main.dataset.00394da4-c9bc-4a9d-9776-db29ca7c0464")
)
after_match_with_GLM_sudo_r <- function(results_Propensity_Matching_cohort_with_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_with_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM+ GFR  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm - only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

Vars <- c("insample glm - only ACE/ARB", "insample glm all variables ")
Pseudo_R_Squared <- c(1 - glm_treat1$deviance / glm_treat1$null.deviance, 1 - glm_treat2$deviance / glm_treat2$null.deviance)

df <- data.frame(Vars, Pseudo_R_Squared)
return(df)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.6b999ec1-5859-4610-86de-e5d5e0441c69"),
    results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.dc40e44c-80a4-4e85-86bf-07b4a9ed206a")
)
after_match_with_LosticModel_sudo_r <- function(results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM+ GFR  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm - only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

Vars <- c("insample glm - only ACE/ARB", "insample glm all variables ")
Pseudo_R_Squared <- c(1 - glm_treat1$deviance / glm_treat1$null.deviance, 1 - glm_treat2$deviance / glm_treat2$null.deviance)

df <- data.frame(Vars, Pseudo_R_Squared)
return(df)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.6662a5a8-3346-4404-af9b-9174e6d86b0e"),
    results_Propensity_Matching_cohort_without_GFR_GLM=Input(rid="ri.foundry.main.dataset.7bbd5171-673f-4bf0-8522-a9f3419123f8")
)
after_match_without_GFR_GLM_ACE_ARB_only <- function(results_Propensity_Matching_cohort_without_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_GLM

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.843f92c9-1c8c-4fe0-8a76-4329f08b9163"),
    results_Propensity_Matching_cohort_without_GFR_GLM=Input(rid="ri.foundry.main.dataset.7bbd5171-673f-4bf0-8522-a9f3419123f8")
)
after_match_without_GFR_GLM_ALL <- function(results_Propensity_Matching_cohort_without_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.64f65335-1f0c-43d7-a32a-cd52d4fb8fcb"),
    results_Propensity_Matching_cohort_without_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.77db6ef3-28d6-4ce8-8ef4-536367d5028d")
)
after_match_without_GFR_LogisticModel_ACE_ARB_only <- function(results_Propensity_Matching_cohort_without_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_LogisticModel

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.0ed5e308-eb26-4241-9aef-b00e092c9ff3"),
    results_Propensity_Matching_cohort_without_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.77db6ef3-28d6-4ce8-8ef4-536367d5028d")
)
after_match_without_GFR_LogisticModel_ALL <- function(results_Propensity_Matching_cohort_without_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_LogisticModel

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

# momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov2[i])
#     print(diff_means_after[i])
# }
# print(diff_means_after[1])

# tab_after <- map_df(diff_means_after, tidy)
# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('ACE_ARB'), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# return(tab)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# tab <- tidy(glm_treat1)
# tab <- tab %>% 
#   rename(
#     # mean_difference = estimate,
#     z_statistic = statistic
#     # df_parameter = parameter,
#     # mean_not_ACE_ARB = estimate1,
#     # Pr_z = p_value
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.3a5241ec-2b96-4894-8870-e990c7b13ce5"),
    results_Propensity_Matching_cohort_without_GFR_GLM=Input(rid="ri.foundry.main.dataset.7bbd5171-673f-4bf0-8522-a9f3419123f8")
)
after_match_without_GLM_sudo_r <- function(results_Propensity_Matching_cohort_without_GFR_GLM) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_GLM

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm - only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

Vars <- c("insample glm - only ACE/ARB", "insample glm all variables ")
Pseudo_R_Squared <- c(1 - glm_treat1$deviance / glm_treat1$null.deviance, 1 - glm_treat2$deviance / glm_treat2$null.deviance)

df <- data.frame(Vars, Pseudo_R_Squared)
return(df)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.dedab404-0424-40e7-ae11-85f556329d8b"),
    results_Propensity_Matching_cohort_without_GFR_LogisticModel=Input(rid="ri.foundry.main.dataset.77db6ef3-28d6-4ce8-8ef4-536367d5028d")
)
after_match_without_LogisticModel_sudo_r <- function(results_Propensity_Matching_cohort_without_GFR_LogisticModel) {

library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

dta_m <- results_Propensity_Matching_cohort_without_GFR_LogisticModel

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))

momeCTR_cov2 <- c('ACE_ARB') # 'GFR' ,
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, 'mort_hospice'] ~dta_m[, v])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov2[i])
    print(diff_means_after[i])
}
print(diff_means_after[1])

tab_after <- map_df(diff_means_after, tidy)
tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('ACE_ARB'), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
print(summary(glm_treat1))

tab <- tidy(glm_treat1)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
print(summary(glm_treat2))

tab <- tidy(glm_treat2)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)

# # # Pseudo R Squared
print("*** Pseudo R Squared")
cat("insample glm - only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

Vars <- c("insample glm - only ACE/ARB", "insample glm all variables ")
Pseudo_R_Squared <- c(1 - glm_treat1$deviance / glm_treat1$null.deviance, 1 - glm_treat2$deviance / glm_treat2$null.deviance)

df <- data.frame(Vars, Pseudo_R_Squared)
return(df)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.4aa101b8-2164-45df-8c5e-42840ba8282c"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_Propensity_Matching_cohort_negative_with_GFR <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

tab_after <- map_df(diff_means_after, tidy)

tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR',  'distance' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
   fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
#    fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 4, widths = c(1, 0.85)
)
# return() ##################
return(tab)

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.dc40e44c-80a4-4e85-86bf-07b4a9ed206a"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_Propensity_Matching_cohort_negative_with_GFR_LogisticModel <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))
return(dta_m)
# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# tab_after <- map_df(diff_means_after, tidy)

# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR',  'distance' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# Checking balance
library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# return() ##################
# return(tab)

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2c4cc79d-2cf7-4187-8c33-eda69cd8eb4d"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_Propensity_Matching_cohort_negative_without_GFR <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

library(broom)
library(purrr)

tab_after <- map_df(diff_means_after, tidy)

tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','distance' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
#    fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 3, widths = c(1, 0.85)
)
# return() ##################
return(tab)

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.ec09b8e6-5e55-49ce-939d-c70d6e7c9dd7"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_Propensity_Matching_cohort_with_GFR <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

tab_after <- map_df(diff_means_after, tidy)

tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR',  'distance' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
   fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
#    fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 4, widths = c(1, 0.85)
)
# return() ##################
return(tab)

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.00394da4-c9bc-4a9d-9776-db29ca7c0464"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_Propensity_Matching_cohort_with_GFR_GLM <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))
return(dta_m)
# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# tab_after <- map_df(diff_means_after, tidy)

# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR',  'distance' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# Checking balance
library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# return() ##################
# return(tab)

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.db6268fd-eef4-447d-bc97-361999e1f62e"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_Propensity_Matching_cohort_without_GFR <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))

# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

library(broom)
library(purrr)

tab_after <- map_df(diff_means_after, tidy)

tab <- as_tibble(tab_after, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','distance' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
print("*** plot propensity vs. mean of each variable")
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$ACE_ARB <- as.factor(dta$ACE_ARB)
  ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw()
}
# Checking balance
library(gridExtra)
grid.arrange(
   fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
   fn_bal(dta_m, "gender"),
   fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
   fn_bal(dta_m, "hx_CHF"),
   fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
   fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
#    fn_bal(dta_m, "GFR_imputed"),
   
   nrow = 3, widths = c(1, 0.85)
)
# return() ##################
return(tab)

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7bbd5171-673f-4bf0-8522-a9f3419123f8"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_Propensity_Matching_cohort_without_GFR_GLM <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))
return(dta_m)
# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

library(broom)
library(purrr)

# tab_after <- map_df(diff_means_after, tidy)

# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','distance' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# Checking balance
library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# return() ##################
# return(tab)

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.77db6ef3-28d6-4ce8-8ef4-536367d5028d"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_Propensity_Matching_cohort_without_GFR_LogisticModel <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# return() #####################################

# omit all incomplete cases
momeCTR_nomiss <- momeCTR %>%  
  select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
  na.omit() 

# Propensity Matching #  GFR  GFR_imputed ADDED
print("*** Propensity Matched")
mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
                     method = "nearest", data = momeCTR_nomiss)
dta_m <- match.data(mod_match)
print(dim(dta_m))
print(summary(mod_match))
return(dta_m)
# Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
for(i in 1:length(diff_means_after)) {
    print(momeCTR_cov[i])
    print(diff_means_after[i])
}
print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

library(broom)
library(purrr)

# tab_after <- map_df(diff_means_after, tidy)

# tab <- as_tibble(tab_after, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','distance' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# Checking balance
library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# return() ##################
# return(tab)

# Estimating treatment effects
print("***  Estimating treatment effects")
trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.90307f2c-9f5e-413f-b202-043b7f595adc"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_before_propensity_Cohort_Negative_with_GFR_LogisticModel_ <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1


momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Linear Model before propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model before propensity matching")
m_ps <- glm( mort_hospice~ ACE_ARB+race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)




tab <- tidy(m_ps)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f0add052-1b7e-4a03-b6b1-48a4bfe7911b"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_before_propensity_Cohort_Negative_without_GFR_LogisticModel_ <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 0,]

# Linear Model before propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model before propensity matching")
m_ps <- glm( mort_hospice~ACE_ARB + race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

tab <- tidy(m_ps)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2c25e12d-5f2e-452a-8b2b-1028d33449cc"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_before_propensity_with_GFR_LogisticModel_ <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Linear Model before propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model before propensity matching")
m_ps <- glm( mort_hospice~ ACE_ARB+race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

tab <- tidy(m_ps)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.4ee41202-b0d3-40de-b71a-662734a47966"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_before_propensity_without_GFR_LogisticModel_ <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 0,]

# Linear Model before propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model before propensity matching")
m_ps <- glm( mort_hospice~ACE_ARB + race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

tab <- tidy(m_ps)
tab <- tab %>% 
  rename(
    # mean_difference = estimate,
    z_statistic = statistic
    # df_parameter = parameter,
    # mean_not_ACE_ARB = estimate1,
    # Pr_z = p_value
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
return(tab)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f423db10-04bd-4296-b52d-48525c5cfd35"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_cohort_negative_with_GFR_ttest <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# library(broom)
# library(purrr)

tab_before <- map_df(diff_means_before, tidy)
# print(tab_before)

tab <- as_tibble(tab_before, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
print("return")
# return (tab)

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

return (tab)

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.75ea261b-0b04-452e-9e25-cc640c704cf4"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_cohort_negative_without_GFR <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)
is.num <- sapply(momeTBL, is.numeric)
momeTBL[is.num] <- lapply(momeTBL[is.num], round, 4)
return(momeTBL)

# # unadjusted difference in mortality and hospice before matching
# print("*** Unadjusted difference in mortality and hospice before matching")
# nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
# print(nm_difference)

# # Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
# print("*** Difference-in-means: pre-treatment covariates")
# momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

# pre_treatment_covariates <- momeCTR %>%
#   group_by(ACE_ARB) %>%
#   select(one_of(momeCTR_cov)) %>%
#   summarise_all(funs(mean(., na.rm = T)))
# print(pre_treatment_covariates)

# # Statistical significance in difference of the means for variable before treatments
# print("*** Statistical significance in difference of the means for variable before treatments")
# diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_before)) {
#     print(momeCTR_cov[i])
#     print(diff_means_before[i])
# }

# library(broom)
# library(purrr)

# tab_before <- map_df(diff_means_before, tidy)
# print(tab_before)

# tab <- as_tibble(tab_before, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# # Linear Model for propensity matching # GFR GFR_imputed ADDED 
# print("*** Linear Model for propensity matching")
# m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
#             family = binomial(), data = momeCTR)
# print(summary(m_ps))
# print(m_ps)

# # output <- new.output()
# # output_fs <- output$fileSystem()
# # pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# # ggsave("Probability_of_prescribed.jpeg")
# # return() ########################################

# # Dataframe of propensity and treatment status
# prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
#                      ACE_ARB = m_ps$model$ACE_ARB)
# print(head(prs_df))

# # histogram of Estimated propensity by treatment
# labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
# hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# # return() #####################################

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.19b0bec4-e31f-4f6d-8a39-e2ff16cbc80d"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_cohort_negative_without_GFR_ttest  <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

library(broom)
library(purrr)

tab_before <- map_df(diff_means_before, tidy)
print(tab_before)

tab <- as_tibble(tab_before, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

# return(tab)

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################
return(tab)
# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.391578cf-67f4-40e5-a947-ff7e89856b4f"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
results_cohort_negetive_with_GFR <- function(ACE_ARB_Pull_negative) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull_negative$gender[ACE_ARB_Pull_negative$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull_negative %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>% 
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))

is.num <- sapply(momeTBL, is.numeric)
momeTBL[is.num] <- lapply(momeTBL[is.num], round, 4)
return(momeTBL)

# # unadjusted difference in mortality and hospice before matching
# print("*** Unadjusted difference in mortality and hospice before matching")
# nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
# print(nm_difference)

# # Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
# print("*** Difference-in-means: pre-treatment covariates")
# momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

# pre_treatment_covariates <- momeCTR %>%
#   group_by(ACE_ARB) %>%
#   select(one_of(momeCTR_cov)) %>%
#   summarise_all(funs(mean(., na.rm = T)))
# print(pre_treatment_covariates)

# # Statistical significance in difference of the means for variable before treatments
# print("*** Statistical significance in difference of the means for variable before treatments")
# diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_before)) {
#     print(momeCTR_cov[i])
#     print(diff_means_before[i])
#     # append(before_treat,)
# }

# # library(broom)
# # library(purrr)

# tab_before <- map_df(diff_means_before, tidy)
# # print(tab_before)

# tab <- as_tibble(tab_before, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','mean_difference'))
# print(tab)
# print("return")
# return (tab)

# # Linear Model for propensity matching # GFR GFR_imputed ADDED 
# print("*** Linear Model for propensity matching")
# m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
#             family = binomial(), data = momeCTR)
# print(summary(m_ps))
# print(m_ps)

# # output <- new.output()
# # output_fs <- output$fileSystem()
# # pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# # ggsave("Probability_of_prescribed.jpeg")
# # return() ########################################

# # Dataframe of propensity and treatment status
# prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
#                      ACE_ARB = m_ps$model$ACE_ARB)
# print(head(prs_df))

# # histogram of Estimated propensity by treatment
# labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
# hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# # return() #####################################

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.d968e9c5-3102-4bda-8963-3ebd18abc828"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_overal <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
#momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients))#,
    # sum_GFR_exist = ( sum(GFR_exist)),
    # sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ) # ,'GFR'

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

tab_before <- map_df(diff_means_before, tidy)
# print(tab_before)

tab <- as_tibble(tab_before, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
print("return")
# return (tab)

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  , # + GFR 
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# # return() #####################################
return (tab)

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   , # + GFR
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance') # 'GFR' ,
# diff_means <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means)) {
#     print(momeCTR_cov[i])
#     print(diff_means[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    #fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.1f03b690-dd89-47c4-9bfb-f6cb4729604d"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_overal_Summary <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
#momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients))#,
    # sum_GFR_exist = ( sum(GFR_exist)),
    # sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)
is.num <- sapply(momeTBL, is.numeric)
momeTBL[is.num] <- lapply(momeTBL[is.num], round, 4)
return(momeTBL)

# # unadjusted difference in mortality and hospice before matching
# print("*** Unadjusted difference in mortality and hospice before matching")
# nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
# print(nm_difference)

# # Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
# print("*** Difference-in-means: pre-treatment covariates")
# momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ) # ,'GFR'

# pre_treatment_covariates <- momeCTR %>%
#   group_by(ACE_ARB) %>%
#   select(one_of(momeCTR_cov)) %>%
#   summarise_all(funs(mean(., na.rm = T)))
# print(pre_treatment_covariates)

# # Statistical significance in difference of the means for variable before treatments
# print("*** Statistical significance in difference of the means for variable before treatments")
# diff_means <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
# for(i in 1:length(diff_means)) {
#     print(momeCTR_cov[i])
#     print(diff_means[i])
# }

# # Linear Model for propensity matching # GFR GFR_imputed ADDED 
# print("*** Linear Model for propensity matching")
# m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM  , # + GFR 
#             family = binomial(), data = momeCTR)
# print(summary(m_ps))
# print(m_ps)

# # output <- new.output()
# # output_fs <- output$fileSystem()
# # pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# # ggsave("Probability_of_prescribed.jpeg")
# # return() ########################################

# # Dataframe of propensity and treatment status
# prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
#                      ACE_ARB = m_ps$model$ACE_ARB)
# print(head(prs_df))

# # histogram of Estimated propensity by treatment
# labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
# hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# # return() #####################################

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   , # + GFR
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance') # 'GFR' ,
# diff_means <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means)) {
#     print(momeCTR_cov[i])
#     print(diff_means[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    #fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m) # + GFR 
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.1ac776a0-9fdf-4440-aa53-94187631cabe"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_with_GFR <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
    # append(before_treat,)
}

# library(broom)
# library(purrr)

tab_before <- map_df(diff_means_before, tidy)
# print(tab_before)

tab <- as_tibble(tab_before, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)
print("return")
# return (tab)

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################

return (tab)

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c91906b2-b659-4315-9237-fc2c0eb448e1"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_with_GFR_summary <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[momeCTR$GFR_exist == 1,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))

is.num <- sapply(momeTBL, is.numeric)
momeTBL[is.num] <- lapply(momeTBL[is.num], round, 4)
return(momeTBL)

# # unadjusted difference in mortality and hospice before matching
# print("*** Unadjusted difference in mortality and hospice before matching")
# nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
# print(nm_difference)

# # Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
# print("*** Difference-in-means: pre-treatment covariates")
# momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' )

# pre_treatment_covariates <- momeCTR %>%
#   group_by(ACE_ARB) %>%
#   select(one_of(momeCTR_cov)) %>%
#   summarise_all(funs(mean(., na.rm = T)))
# print(pre_treatment_covariates)

# # Statistical significance in difference of the means for variable before treatments
# print("*** Statistical significance in difference of the means for variable before treatments")
# diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_before)) {
#     print(momeCTR_cov[i])
#     print(diff_means_before[i])
#     # append(before_treat,)
# }

# # library(broom)
# # library(purrr)

# tab_before <- map_df(diff_means_before, tidy)
# # print(tab_before)

# tab <- as_tibble(tab_before, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM','GFR' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','mean_difference'))
# print(tab)
# print("return")
# return (tab)

# # Linear Model for propensity matching # GFR GFR_imputed ADDED 
# print("*** Linear Model for propensity matching")
# m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM + GFR  ,
#             family = binomial(), data = momeCTR)
# print(summary(m_ps))
# print(m_ps)

# # output <- new.output()
# # output_fs <- output$fileSystem()
# # pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# # ggsave("Probability_of_prescribed.jpeg")
# # return() ########################################

# # Dataframe of propensity and treatment status
# prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
#                      ACE_ARB = m_ps$model$ACE_ARB)
# print(head(prs_df))

# # histogram of Estimated propensity by treatment
# labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
# hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# # return() #####################################

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM', 'GFR' ,  'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
#    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 4, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM + GFR  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7c6d630d-5685-4610-9ce5-d7c94d166d73"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_without_GFR <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)

# unadjusted difference in mortality and hospice before matching
print("*** Unadjusted difference in mortality and hospice before matching")
nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
print(nm_difference)

# Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
print("*** Difference-in-means: pre-treatment covariates")
momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

pre_treatment_covariates <- momeCTR %>%
  group_by(ACE_ARB) %>%
  select(one_of(momeCTR_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))
print(pre_treatment_covariates)

# Statistical significance in difference of the means for variable before treatments
print("*** Statistical significance in difference of the means for variable before treatments")
diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
for(i in 1:length(diff_means_before)) {
    print(momeCTR_cov[i])
    print(diff_means_before[i])
}

library(broom)
library(purrr)

tab_before <- map_df(diff_means_before, tidy)
print(tab_before)

tab <- as_tibble(tab_before, col_labels = TRUE)
tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ), tab)
tab <- tab %>% 
  rename(
    mean_difference = estimate,
    t_statistic = statistic,
    df_parameter = parameter,
    mean_not_ACE_ARB = estimate1,
    mean_ACE_ARB = estimate2
    )
is.num <- sapply(tab, is.numeric)
tab[is.num] <- lapply(tab[is.num], round, 8)
tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
print(tab)

# return(tab)

# Linear Model for propensity matching # GFR GFR_imputed ADDED 
print("*** Linear Model for propensity matching")
m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
            family = binomial(), data = momeCTR)
print(summary(m_ps))
print(m_ps)

# output <- new.output()
# output_fs <- output$fileSystem()
# pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
plot(predict(m_ps),residuals(m_ps))
# ggsave("Probability_of_prescribed.jpeg")
# return() ########################################

# Dataframe of propensity and treatment status
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     ACE_ARB = m_ps$model$ACE_ARB)
print(head(prs_df))

# histogram of Estimated propensity by treatment
labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

gplt <- ggplot(hist, aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~ACE_ARB) +
    xlab("Probability of being prescribed medicine of interest") +
    theme_bw()
plot(gplt)
# return() #####################################
return(tab)
# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.a6d56517-fb9d-4ca8-96ff-8fc5924fadfc"),
    ACE_ARB_Pull=Input(rid="ri.foundry.main.dataset.a038d4e4-8900-4748-b40d-2daef98beeb7")
)
results_cohort_without_GFR_summary <- function(ACE_ARB_Pull) {

# R packages required
library(MatchIt)
library(dplyr)
library(ggplot2)

# Hamid # ROOT cause of problem for lapply # gender was Null (NA)
# print(head(ACE_ARB_Pull %>% select('gender','gender_pat')))
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "FEMALE"] <- 0
ACE_ARB_Pull$gender[ACE_ARB_Pull$gender_pat == "MALE"] <- 1

momeCTR  <- ACE_ARB_Pull %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN, GFR) %>%
mutate(x = as.numeric(ACE_ARB))

# imputation
momeCTR <- momeCTR %>% mutate(GFR_exist = ifelse(is.na(GFR), 0, 1))
# momeCTR$GFR[is.na(momeCTR$GFR)]<-mean(momeCTR$GFR,na.rm=TRUE)
print(head(momeCTR))
momeCTR <- momeCTR[ momeCTR$GFR_exist == 0,]

# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Dataset")
summ <- momeCTR %>%
    summarise(n_patients = n(),
    sum_morth = sum(mort_hospice),
    std_error = (sd(mort_hospice) / sqrt(n_patients)),
    sum_GFR_exist = ( sum(GFR_exist)),
    sum_GFR_missing = n() -sum(GFR_exist)
    )
print(summ)

momeCTR  <- momeCTR  %>% dplyr::select(visit_concept_name, Severity_Type, age_pat, ACE_ARB, race_AA, Ethnicity, gender, mort_hospice, hx_DM, hx_CHF, hx_HTN)
# Pre-analysis using non-matched data of mortality and hospice outcome
print("*** Pre-analysis using non-matched data of mortality and hospice outcome")
momeTBL <- momeCTR %>%
  group_by(ACE_ARB) %>%
  summarise(n_patients = n(),
            mean_morth = mean(mort_hospice),
            std_error = sd(mort_hospice) / sqrt(n_patients))
print(momeTBL)
is.num <- sapply(momeTBL, is.numeric)
momeTBL[is.num] <- lapply(momeTBL[is.num], round, 4)
return(momeTBL)

# # unadjusted difference in mortality and hospice before matching
# print("*** Unadjusted difference in mortality and hospice before matching")
# nm_difference <- with(momeCTR, t.test(mort_hospice ~ ACE_ARB))
# print(nm_difference)

# # Difference-in-means: pre-treatment covariates # GFR GFR_imputed ADDED
# print("*** Difference-in-means: pre-treatment covariates")
# momeCTR_cov <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' )

# pre_treatment_covariates <- momeCTR %>%
#   group_by(ACE_ARB) %>%
#   select(one_of(momeCTR_cov)) %>%
#   summarise_all(funs(mean(., na.rm = T)))
# print(pre_treatment_covariates)

# # Statistical significance in difference of the means for variable before treatments
# print("*** Statistical significance in difference of the means for variable before treatments")
# diff_means_before <- lapply(momeCTR_cov, function(v){t.test(momeCTR[, v] ~ momeCTR[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_before)) {
#     print(momeCTR_cov[i])
#     print(diff_means_before[i])
# }

# library(broom)
# library(purrr)

# tab_before <- map_df(diff_means_before, tidy)
# print(tab_before)

# tab <- as_tibble(tab_before, col_labels = TRUE)
# tab <- cbind(covariate = c('race_AA', 'age', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM' ), tab)
# tab <- tab %>% 
#   rename(
#     mean_difference = estimate,
#     t_statistic = statistic,
#     df_parameter = parameter,
#     mean_not_ACE_ARB = estimate1,
#     mean_ACE_ARB = estimate2
#     )
# is.num <- sapply(tab, is.numeric)
# tab[is.num] <- lapply(tab[is.num], round, 8)
# tab <- tab %>% select(,-c('method','alternative','t_statistic','df_parameter')) # ,'mean_difference'
# print(tab)
# print("return")
# return(tab)

# # Linear Model for propensity matching # GFR GFR_imputed ADDED 
# print("*** Linear Model for propensity matching")
# m_ps <- glm(ACE_ARB ~ race_AA + age_pat + gender + hx_CHF + hx_HTN + hx_DM ,
#             family = binomial(), data = momeCTR)
# print(summary(m_ps))
# print(m_ps)

# # output <- new.output()
# # output_fs <- output$fileSystem()
# # pdf(output_fs$get_path("Probability_of_prescribed.pdf", 'w'))
# plot(predict(m_ps),residuals(m_ps))
# # ggsave("Probability_of_prescribed.jpeg")
# # return() ########################################

# # Dataframe of propensity and treatment status
# prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
#                      ACE_ARB = m_ps$model$ACE_ARB)
# print(head(prs_df))

# # histogram of Estimated propensity by treatment
# labs <- paste("Actual Medicine Prescribed:", c("Medicine", "Not_Medicine"))
# hist <- prs_df %>% mutate(ACE_ARB = ifelse(ACE_ARB == 1, labs[1], labs[2])) #%>%

# gplt <- ggplot(hist, aes(x = pr_score)) +
#     geom_histogram(color = "white") +
#     facet_wrap(~ACE_ARB) +
#     xlab("Probability of being prescribed medicine of interest") +
#     theme_bw()
# plot(gplt)
# # return() #####################################

# # omit all incomplete cases
# momeCTR_nomiss <- momeCTR %>%  
#   select(mort_hospice, ACE_ARB, one_of(momeCTR_cov)) %>%
#   na.omit() 

# # Propensity Matching #  GFR  GFR_imputed ADDED
# print("*** Propensity Matched")
# mod_match <- matchit(ACE_ARB ~ race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM   ,
#                      method = "nearest", data = momeCTR_nomiss)
# dta_m <- match.data(mod_match)
# print(dim(dta_m))
# print(summary(mod_match))

# # Statistical significance in difference of the means for variable >>>after<<< adjustment # GFR  GFR_imputed ADDED
# print("*** Statistical significance in difference of the means for variable >>>AFTER<<< adjustment")
# momeCTR_cov2 <- c('race_AA', 'age_pat', 'gender', 'hx_HTN', 'hx_CHF', 'hx_DM',   'distance')
# diff_means_after <- lapply(momeCTR_cov2 , function(v){t.test(dta_m[, v] ~dta_m[, 'ACE_ARB'])})
# for(i in 1:length(diff_means_after)) {
#     print(momeCTR_cov[i])
#     print(diff_means_after[i])
# }
# print(dta_m %>% group_by(ACE_ARB) %>% summarise_all(.funs = c(mean="mean")))

# # plot propensity vs. mean of each variable
# print("*** plot propensity vs. mean of each variable")
# fn_bal <- function(dta, variable) {
#   dta$variable <- dta[, variable]
#   dta$ACE_ARB <- as.factor(dta$ACE_ARB)
#   ggplot(dta, aes(x = distance, y = variable, color = ACE_ARB)) +
#     geom_point(alpha = 0.2, size = 1.5) +
#     geom_smooth(method = "loess", se = F) +
#     xlab("Propensity score") +
#     ylab(variable) +
#     theme_bw()
# }
# # Checking balance
# library(gridExtra)
# grid.arrange(
#    fn_bal(dta_m, "age_pat") + theme(legend.position = "none"),
#    fn_bal(dta_m, "gender"),
#    fn_bal(dta_m, "hx_HTN") + theme(legend.position = "none"),
#    fn_bal(dta_m, "hx_CHF"),
#    fn_bal(dta_m, "hx_DM") + theme(legend.position = "none"),
#    fn_bal(dta_m, "race_AA"),
# #    fn_bal(dta_m, "GFR") + theme(legend.position = "none"),
# #    fn_bal(dta_m, "GFR_imputed"),
   
#    nrow = 3, widths = c(1, 0.85)
# )
# # return() ##################

# # Estimating treatment effects
# print("***  Estimating treatment effects")
# trt_t <- with(dta_m, t.test(mort_hospice ~ ACE_ARB))
# print(trt_t)

# glm_treat1 <-glm(formula = mort_hospice ~ ACE_ARB, family = binomial, data = dta_m)
# print(summary(glm_treat1))

# glm_treat2 <- glm(mort_hospice ~ ACE_ARB + race_AA + age_pat + gender + hx_HTN + hx_CHF + hx_DM  , family = binomial, data = dta_m)
# print(summary(glm_treat2))

# # plotting Residuals
# print("*** plotting Residuals")
# plot(predict(glm_treat2),residuals(glm_treat2))

# # Pseudo R Squared
# print("*** Pseudo R Squared")
# cat("insample glm only ACE/ARB ",1 - glm_treat1$deviance / glm_treat1$null.deviance)
# cat("\ninsample glm all variables ",1 - glm_treat2$deviance / glm_treat2$null.deviance)

}

#

@transform_pandas(
    Output(rid="ri.vector.main.execute.628260ba-15bd-4843-91ea-a3dade0e29c2"),
    ACE_ARB_Pull_negative=Input(rid="ri.foundry.main.dataset.7ac258ca-9aeb-40d4-9f17-3a9b1de485e4")
)
unnamed_4 <- function(ACE_ARB_Pull_negative) {
    
}


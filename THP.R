---
title: "THP Global Priorities - Presentation Results"
author: "Cedric Antunes"
date: "September 2025"
output: pdf_document
---

```{r setup, include=FALSE}
rm(list = ls())

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(echo = FALSE)

library(dplyr)
library(tidyr)
library(fixest)
library(vtable)
library(ggthemes)
library(purrr)
library(gt)
library(stringr)
library(ggplot2)
library(anthro)
```

```{r}
#### functions #### 
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  
  results <- purrr::map_dfr(outcomes, function(outcome) {
    
    
    if (is.null(controls)) {
      fml <- as.formula(paste0(outcome, " ~ treat"))
    } else {
      fml <- as.formula(paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    }
    
    
    model <- tryCatch({
      feols(fml, data = data, se = "hetero", weights = weights)
    }, error = function(e) return(NULL))
    
    
    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]
      se_val <- se(model)["treat"]
      
      tibble(
        outcome = outcome,
        estimate = est,
        std_error = se_val,
        ci_low = est - 1.96 * se_val,
        ci_high = est + 1.96 * se_val,
        est_type = model_type, 
      )
    } else {
      tibble(
        outcome = outcome,
        estimate = NA,
        std_error = NA,
        ci_low = NA,
        ci_high = NA,
        est_type = model_type
      )
    }
  })
  
  return(results)
}

plot_by_section <- function(section_name, results_df) {
  plot_data <- results_df %>% filter(section == section_name)
  
  ggplot(plot_data, aes(x = reorder(label_wrapped, estimate), y = estimate, color = est_type, shape = est_type)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, position = position_dodge(width = 0.6)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "Outcome", y = "Treatment Effect", title = paste("Regression Results:", section_name)) +
    scale_color_manual(values = c("Bivariate" = "#F28968", "Controls" = "#6CBF84", "Controls  (+Weights)" = '#A0522D')) +
    scale_shape_manual(values = c("Bivariate" = 1, "Controls" = 2, "Controls  (+Weights)" = 3)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom") + 
    labs(shape = "Estimation Type") + 
    labs(color = "Estimation Type")
}
```

```{r, echo = F, warning=FALSE, message=FALSE}
#### Mexico #### 

thp_mex = haven::read_spss('F:/cedric.antunes/Documents/Evaluasi/THP_data_input/DATAOPM_Mexico_database.sav') %>% 
  filter(q16 == 1 | q19 >0) %>% 
  mutate_if(is.numeric, ~ if_else(. %in% c(-99, -98), NA,  .)) %>% 
  mutate_if(is.numeric, ~ if_else(. %in% c(99, 98), NA,  .)) %>% 
  mutate(treat = if_else(e3 == 1, 1,0)) %>% 
  mutate(across(c(q2, q12, q13), ~ as.factor(.))) %>% 
  mutate(q20 = if_else(q20 %in% c(1, 2), 1,0 )) %>% 
  mutate(
    q5 = if_else(q5 %in% c('96', '99'), NA, as.numeric(q5))) %>% 
  mutate(roof = if_else(q9 > 3, 1,0),
         floor = if_else(q8  > 4, 1,0), 
         water = if_else(q10  %in% c(6, 3, 5), 1,0),
         phone = if_else(q11 %in% c(2, 3, 4, 5), 1,0)
  ) %>% 
  mutate(wlth = prcomp(select(., roof, floor, water, phone))$x[,1]) %>% 
  mutate(q18_w = as.numeric(q18)/4) %>% 
  mutate(q20 = if_else(q20 %in% c('1', '2'), 1, 0),
         q28 = if_else(q28 %in% c(96, 99, 100), NA, q28),
         correct_q32 = if_else(q32 == 3, 1, 0),
         correct_q33 = if_else(q33 == 6, 1, 0), 
         breastfeed_knowledge = (correct_q32 + correct_q33) / 2,
         wash_index = (q68 + q69) / 12, 
         norms_index = (q64+q65+(2-q66)/5)
  ) %>% 
  mutate(across(c(q2, q12, q13), ~ as.factor(.))) %>% 
  mutate(q20 = if_else(q20 %in% c(1, 2), 1,0 )) %>% 
  mutate(
    q5 = if_else(q5 %in% c('96', '99'), NA, as.numeric(q5))) %>% 
  mutate(roof = if_else(q9 > 3, 1,0),
         floor = if_else(q8  > 4, 1,0), 
         water = if_else(q10  %in% c(6, 3, 5), 1,0),
         phone = if_else(q11 %in% c(2, 3, 4, 5), 1,0)
  ) %>% 
  mutate(wlth = prcomp(select(., roof, floor, water, phone))$x[,1]) %>% 
  mutate(q24 = q24_1/24 + q24_2)

thp_mex$q45wt <- 2 * rowSums(thp_mex[, c("q45_3", "q45_4", "q45_1", "q45_2")], na.rm = TRUE) +
  1 * rowSums(thp_mex[, c("q45_5", "q45_6")], na.rm = TRUE) +
  4 * rowSums(thp_mex[, c("q45_7", "q45_8", "q45_9")], na.rm = TRUE) +
  0.5 * rowSums(thp_mex[, c("q45_12", "q45_13", "q45_14", "q45_15", 
                            "q45_16",
                            "q45_17", "q45_18", "q45_19", 
                            "q45_20", "q45_21")], na.rm = TRUE)

thp_mex$q48wt <- 2 * rowSums(thp_mex[, c("q48b_3", "q48b_4", "q48b_1", "q48b_2")], na.rm = TRUE) +
  1 * rowSums(thp_mex[, c("q48b_5", "q48b_6")], na.rm = TRUE) +
  4 * rowSums(thp_mex[, c("q48b_7", "q48b_8", "q48b_9")], na.rm = TRUE) +
  0.5 * rowSums(thp_mex[, c("q48b_12", "q48b_13", "q48b_14", "q48b_15", 
                            "q48b_16",
                            "q48b_17", "q48b_18", "q48b_19", 
                            "q48b_20", "q48b_21", "q48b_22", "q48b_23")], na.rm = TRUE)


control_vars = c('q2', 'q3',  'q6', 
                 'q7', 'q12',  'wlth')
```

```{r, echo = F, warning=FALSE, message=FALSE}
vtable::st(thp_mex %>%  select(treat, control_vars), group = 'treat', group.test = TRUE, 
           title = 'Covariate Balance - Full Sample')
```

```{r, echo = F, warning=FALSE, message=FALSE}

fml <- as.formula(paste0('treat ~',  paste(control_vars, collapse = " + ")))

pred_treat = glm(fml, data = thp_mex, 
                 family = binomial(link = 'logit'))

thp_mex$prob = predict(pred_treat, type = 'response')
thp_mex = thp_mex %>%  
  filter(prob > .1 & prob < .9)

vtable::st(thp_mex %>%  select(treat, control_vars), group = 'treat', group.test = TRUE, 
            title = 'Covariate Balance - Trimmed')

thp_mex = thp_mex %>%
  mutate(across(everything(), ~ if_else(. == 99, NA, .)))

thp_scale <- thp_mex %>% 
  mutate(across(
    c('q65', 'q66', 'wash_index', 
      'q19', 'q20', 'q24', 'q25', 'q26', 
      'q27', 'q28', 'q30', 'q31', 'breastfeed_knowledge', 
      'q34', 'q36', 'q37', 'q40', 'q42', 
      'q45wt', 'q48wt', 'q46', 'q49',
      'q64', 'q65', 'q66', 'norms_index',
      't18','t19','t20','t21','t22','t23','t24','t25','t26',
      't27','t28','t29','t30','t31','t32','t33','t34'),
    ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
  ))

thp_scale = thp_scale %>%  
    mutate(weight = if_else(treat == 1, mean(treat) / prob, (1-mean(treat)) / (1 - prob)))

thp_scale = thp_scale %>%  
  filter(prob > .1 & prob < .75)
```


```{r, echo = F, warning=FALSE, message=FALSE}
ggplot(thp_scale, aes(x = prob, fill = factor(treat))) +
  geom_histogram(aes(prob, fill = factor(treat)),
                 position = "identity", 
                 alpha = 0.7) + 
  labs(title = "Common Support by Group (E3 Treatment)", 
       x = "Pr(Treat | X)", 
       y = "Frequency", 
       fill = 'Beneficiary') +
  theme_minimal() +
  scale_fill_manual(values = c("1" = "#F28968", "0" = "#6CBF84")) +
  theme(legend.position = 'bottom') + 
  theme_minimal()
```

```{r, echo = F, warning=FALSE, message=FALSE}

outcome_groups <- list(
  Breastfeeding = c('q19', 'q20', 'q24', 'q25', 'q26', 
                    'q27', 'q28', 'q30', 'q31', 'breastfeed_knowledge'),
  
  Pregnancy = c('q34', 'q36', 'q37', 'q40', 'q42'),
  
  Food_Supplements = c('q45wt', 'q48wt', 'q46', 'q49'),
  
  WASH = c('q68', 'q69', 'wash_index'), 
  
  Knowledge_Skills = c('t18', 't19', 't20', 
                       't21', 't22' , 
                       't23', 't24', 't25', 't26', 
                       't27', 't28', 't29', 't30', 
                       't31','t32', 't33', 't34'),
  
  norms = c('q64', 'q65','q66','norms_index')
  
)

variable_labels <- tibble(
  outcome = unlist(outcome_groups),  
  label = unlist(c(
    c(
      "Children Alive",
      "Has Health Card",
      "Time to First Breastfeed",
      "Non-Breastfeeding in First 2 Days",
      "Still Breastfeeding",
      "Breastfed Times in Last 24 Hrs",
      "Age When Supplemented",
      "Planned Months Exclusive BF",
      "Planned Months Supplementary BF",
      "Breastfeeding Knowledge"
    ),
    
    c(
      "Did you ever receive training on good pregnancy habits in last year?",
      "Months Pregnant at First ANC Visit",
      "Number of ANC Visits During Pregnancy",
      "Took Nutritional Supplements During Pregnancy",
      "Days Took Tablets or Syrup During Pregnancy"
    ),
    
    c(
      "Child Food Consumption Score (FCS)",
      "Respondent Food Consumption Score (FCS)",
      "Nutrition Supplements (Adults)",
      "Nutrition Supplements (Children)"
    ),
    
    c(
      "Wash Hands: Cooking",
      "Wash Hands: Eating",
      "Wash Hands: Index"
    ), 
    
    
    c(
      "Breastfeeding is not necessary",
      "Important to complement breastfeeding",
      "Mothers need extra food",
      "A nutritious diet includes candy",
      "Dietary diversity is crucial",
      "Taking MMS nutritional supplements is good",
      "Malnourished children can improve their growth by taking LNS",
      "Children between 6 and 24 months do not require anything other than breastmilk",
      "Children between 0 and 6 months do not require anything other than breastmilk",
      "Immunization (vaccines) for babies are essential",
      "Using latrines much better than open defecation.",
      "Raising poultry and gardening fruits and vegetables only one month per year is enough to make sure household has good nutririon.",
      "I know how to preserve and use the nutritional supplements",
      "Confidence in breastfeeding",
      "Confidence in feeding baby",
      "Confidence in identifying malnourishment",
      "If I detect a community need, I am able to politely raise it with the relevant government authorities"
    ),
    
    c(
      "Agree: Supplements Means Failure as Caregiver",
      "Agree: Supplements Cause Large Baby & Painful Birth",
      "Confidence in Positioning & Attaching Baby for Breastfeeding",
      "Norms & Stigma Index"
    )
  )))

results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls", 
                        weights = thp_scale$weight) %>% mutate(est_type = paste(est_type, " (+Weights)"))
  bind_rows(bivar, covs, covs_wt) %>% mutate(section = section_name)
})

results_all <- results_all %>%
  left_join(variable_labels, by = "outcome") %>%
  mutate(label_wrapped = stringr::str_wrap(label, width = 25))
```

```{r}
# ------------------------------------------------------------------------------
# Helper function based on Don's analysis --------------------------------------
# ------------------------------------------------------------------------------
# 1) RE-RUN this definition
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- se(model)["treat"]
      tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = nobs(model)
      )
    } else {
      tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
             ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
             control_mean = NA_real_, N = NA_integer_)
    }
  })
}

# 2) REBUILD results_all (rerun your imap_dfr block)
results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>% mutate(est_type = paste(est_type, " (+Weights)"))
  bind_rows(bivar, covs, covs_wt) %>% mutate(section = section_name)
}) %>%
  left_join(variable_labels, by = "outcome") %>%
  mutate(label_wrapped = stringr::str_wrap(label, width = 25))

format_thp_table <- function(results_all, type_pattern, add_labels = TRUE){
  df <- results_all %>%
    dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE))

  # guard: add N if absent (shouldn’t happen if you rebuilt results_all)
  if (!"N" %in% names(df)) df <- dplyr::mutate(df, N = NA_integer_)

  df %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 2),
      se        = round(se, 2),
      ci        = paste0("[", round(conf.low, 2), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 2),
      outcome   = dplyr::if_else(!is.na(label), label, outcome)
    ) %>%
    dplyr::select(outcome, att, se, ci, n, ctrl_mean)
}

mex_table_bivar    <- format_thp_table(results_all, type_pattern = "^bivariate$")
mex_table_controls <- format_thp_table(results_all, type_pattern = "^controls$")
mex_table_weighted <- format_thp_table(results_all, type_pattern = "controls\\s+\\(\\+weights\\)")

print(mex_table_bivar, n = Inf, width = Inf)
print(mex_table_controls, n = Inf, width = Inf)
print(mex_table_weighted, n = Inf, width = Inf)

write.csv(mex_table_bivar,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_bivar.csv")
write.csv(mex_table_controls,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_controls.csv")
```

```{r mexico-restandardization}
# --- build final analysis sample first (after your existing steps) -----------
thp_final <- thp_mex |>
  mutate(across(everything(), ~ if_else(. == 99, NA, .))) %>%
  filter(prob > .1, prob < .75) %>%  # <- your final common-support rule
  mutate(
    weight = if_else(
      treat == 1, mean(treat) / prob, (1 - mean(treat)) / (1 - prob)
    )
  )

# --- variables to standardize: exactly what you model/plot -------------------
vars_to_std <- unique(unlist(outcome_groups))

# --- compute control-only moments on the final sample ------------------------
ctrl_means <- thp_final |>
  filter(treat == 0) |>
  summarise(across(all_of(vars_to_std), ~ mean(., na.rm = TRUE))) |>
  as.list() |> 
  unlist()

ctrl_sds <- thp_final |>
  filter(treat == 0) |>
  summarise(across(all_of(vars_to_std), ~ sd(., na.rm = TRUE))) |>
  as.list() |> 
  unlist()

# guard against zero/NA SDs to avoid divide-by-zero; leave those variables centered-only
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

# --- apply control-relative standardization ----------------------------------
# Option A (recommended): Z-score with control mean & control SD
thp_scale <- thp_final |>
  mutate(across(
    all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()],
    .names = "{.col}"
  ))

# -------------------------------------------------------------------
# Run models on the control-centered thp_scale
# -------------------------------------------------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- se(model)["treat"]
      tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = nobs(model)
      )
    } else {
      tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
             ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
             control_mean = NA_real_, N = NA_integer_)
    }
  })
}

# -------------------------------------------------------------------
# Rebuild results_all using thp_scale (control-centered outcomes)
# -------------------------------------------------------------------
results_all <- purrr::imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>% dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>% dplyr::mutate(section = section_name)
}) %>%
  dplyr::left_join(variable_labels, by = "outcome") %>%
  dplyr::mutate(label_wrapped = stringr::str_wrap(label, width = 25))

# -------------------------------------------------------------------
# Format tables (control_mean should be ~0 now for each outcome)
# -------------------------------------------------------------------
format_thp_table <- function(results_all, type_pattern, add_labels = TRUE){
  df <- results_all %>%
    dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE))

  if (!"N" %in% names(df)) df <- dplyr::mutate(df, N = NA_integer_)

  df %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 2),
      se        = round(se, 2),
      ci        = paste0("[", round(conf.low, 2), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 2),
      outcome   = dplyr::if_else(!is.na(label), label, outcome)
      # Note: att is in "control-SD units" because outcomes were standardized on the control group.
    ) %>%
    dplyr::select(outcome, att, se, ci, n, ctrl_mean)
}

mex_table_bivar    <- format_thp_table(results_all, type_pattern = "^bivariate$")
mex_table_controls <- format_thp_table(results_all, type_pattern = "^controls$")
mex_table_weighted <- format_thp_table(results_all, type_pattern = "controls\\s+\\(\\+weights\\)")

print(mex_table_bivar, n = Inf, width = Inf)
print(mex_table_controls, n = Inf, width = Inf)
print(mex_table_weighted, n = Inf, width = Inf)

write.csv(mex_table_bivar,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_bivar_revised.csv")
write.csv(mex_table_controls,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_controls_revised.csv")
write.csv(mex_table_weighted,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_controls_weighted_revised.csv")
```

```{r mexico-all-regressors}
# -------------------------------------------------------------------
# 1) Outcomes: include index components, keep raw var names
# -------------------------------------------------------------------
outcome_groups <- list(
  Breastfeeding = c('q19','q20','q24','q25','q26','q27','q28','q30','q31',
                    'breastfeed_knowledge', 'correct_q32', 'correct_q33'),
  Pregnancy     = c('q34','q36','q37','q40','q42'),
  Food_Supplements = c('q45wt','q48wt','q46','q49'),
  WASH          = c('q68','q69','wash_index'),
  Knowledge_Skills = c('t18','t19','t20','t21','t22','t23','t24','t25','t26',
                       't27','t28','t29','t30','t31','t32','t33','t34'),
  norms         = c('q64','q65','q66','norms_index')
)

# -------------------------------------------------------------------
# 2) Build final sample + control-relative standardization
#    (keeps SPSS variable labels on columns where they exist)
# -------------------------------------------------------------------
thp_final <- thp_mex %>%
  dplyr::mutate(across(everything(), ~ if_else(. == 99, NA, .))) %>%
  dplyr::filter(prob > .1, prob < .75) %>%
  dplyr::mutate(weight = if_else(treat == 1, mean(treat) / prob, (1 - mean(treat)) / (1 - prob)))

vars_to_std <- unique(unlist(outcome_groups))

ctrl_means <- thp_final %>%
  dplyr::filter(treat == 0) %>%
  dplyr::summarise(across(all_of(vars_to_std), ~ mean(., na.rm = TRUE))) %>%
  as.list() |> unlist()

ctrl_sds <- thp_final %>%
  dplyr::filter(treat == 0) %>%
  dplyr::summarise(across(all_of(vars_to_std), ~ sd(., na.rm = TRUE))) %>%
  as.list() |> unlist()
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

thp_scale <- thp_final %>%
  dplyr::mutate(across(
    all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()]
  ))

# -------------------------------------------------------------------
# 3) Run models (unchanged)
# -------------------------------------------------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(fixest::feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- fixest::se(model)["treat"]
      tibble::tibble(
        outcome      = outcome,                  # raw var name
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = stats::nobs(model)
      )
    } else {
      tibble::tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
                     ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
                     control_mean = NA_real_, N = NA_integer_)
    }
  })
}

results_all <- purrr::imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>%
             dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>%
    dplyr::mutate(section = section_name)
})

# -------------------------------------------------------------------
# 4) Format tables: keep raw var names; no relabelling
# -------------------------------------------------------------------
format_thp_table <- function(results_all, type_pattern){
  df <- results_all %>% dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE))
  if (!"N" %in% names(df)) df <- dplyr::mutate(df, N = NA_integer_)
  df %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 2),
      se        = round(se, 2),
      ci        = paste0("[", round(conf.low, 2), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 2)
    ) %>%
    dplyr::select(section, outcome, att, se, ci, n, ctrl_mean)
}

mex_table_bivar_all_regressors    <- format_thp_table(results_all, "^bivariate$")
mex_table_controls_all_regressors <- format_thp_table(results_all, "^controls$")
mex_table_weighted_all_regressors <- format_thp_table(results_all, "controls\\s+\\(\\+weights\\)")

print(mex_table_bivar, n = Inf, width = Inf)
print(mex_table_controls, n = Inf, width = Inf)
print(mex_table_weighted, n = Inf, width = Inf)

# Optional: write CSVs
write.csv(mex_table_bivar_all_regressors,    
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_bivar_all_regressors.csv",             
          row.names = FALSE)

write.csv(mex_table_controls_all_regressors, 
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_controls_all_regressors.csv",          
          row.names = FALSE)

write.csv(mex_table_weighted_all_regressors,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/mex_controls_weights_all_regressors.csv",  
          row.names = FALSE)
```

\newpage 

# THP Bangladesh 

```{r}
thp_raw <- haven::read_sav('F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_Clean_Data.sav')

thp_clean <- thp_raw %>%
  mutate(across(where(is.numeric), ~ if_else(. %in% c(99, 999, 9999), NA_real_, .))) %>%
  mutate(
    
    yes_e3 = if_else(e3 == 1, 1, 0),
    
    
    q20 = if_else(q20 %in% c('1', '2'), 1, 0),
    q24 = if_else(q24 == 999, NA_real_, q24),
    q28 = if_else(q28 %in% c(96, 99, 100), NA, q28),
    
    
    correct_q32 = if_else(q32 == 3, 1, 0),
    correct_q33 = if_else(q33 == 6, 1, 0),
    breastfeed_knowledge = (correct_q32 + correct_q33) / 2,
    
    wash_index = (q65 + q66) / 12
  )

thp_clean$q43wt <- 2 * rowSums(thp_clean[, c("q43_03", "q43_04", "q43_01", "q43_02")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q43_05", "q43_06")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q43_07", "q43_08", "q43_09")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q43_012", "q43_013", "q43_014", "q43_015", "q43_016",
                              "q43_017", "q43_018", "q43_019", "q43_020", "q43_021")], na.rm = TRUE)

thp_clean$q46wt <- 2 * rowSums(thp_clean[, c("q46_03", "q46_04", "q46_01", "q46_02")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q46_05", "q46_06")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q46_07", "q46_08", "q46_09")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q46_012", "q46_013", "q46_014", "q46_015", "q46_016",
                              "q46_017", "q46_018", "q46_019", "q46_020", "q46_021")], na.rm = TRUE)

thp_clean <- thp_clean %>% 
  mutate(across(c(q2, q12, q13), ~ as.factor(.))) %>% 
  mutate(q20 = if_else(q20 %in% c(1, 2), 1,0 )) %>% 
  mutate(
    q5 = if_else(q5 %in% c('96', '99'), NA, as.numeric(q5))) %>% 
  mutate(roof = if_else(q9 > 3, 1,0),
         floor = if_else(q8  > 4, 1,0), 
         water = if_else(q10  %in% c(6, 3, 5), 1,0),
         phone = if_else(q11 %in% c(2, 3, 4, 5), 1,0)
  ) %>% 
  mutate(wlth = prcomp(select(., roof, floor, water, phone))$x[,1]) %>% 
  mutate(q18_w = as.numeric(q18)/4, 
         q21_w = as.numeric(as.Date(b13)-as.Date(q21))/7) %>% 
  mutate(yes_e3 = as.numeric(e3 == 1, 1,0)) 
thp_clean$treat = thp_clean$yes_e3

fml <- as.formula(paste0('yes_e3 ~',  paste(control_vars, collapse = " + ")))
```

```{r, echo = F, warning=FALSE, message=FALSE}

outcome_groups <- list(
  Breastfeeding = c('q19', 'q20', 'q24', 'q25', 'q26', 
                    'q27', 'q28', 'q30', 'q31', 'breastfeed_knowledge'),
  
  Pregnancy = c('q34', 'q36', 'q37', 'q40', 'q42'),
  
  Food_Supplements = c('q43wt','q46wt', 'q44', 'q47'),
  
  WASH = c('q65', 'q66', 'wash_index')
  
)

control_vars = c('q2', 'q3', 'q5', 'q6', 
                 'q7', 'q12', 'q13', 'wlth', 'q16')

thp_scale <- thp_clean %>% 
  mutate(across(
    c('q65','q66','wash_index',
      'q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge',
      'q34','q36','q37','q40','q42',
      'q43wt','q46wt','q44','q47'),
    ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
  ))

thp_scale$treat  <- thp_scale$yes_e3
thp_scale$weight <- ifelse(thp_scale$treat == 1, 1 / thp_scale$prob, 1 / (1 - thp_scale$prob))

results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls", 
                        weights = thp_scale$weight) %>% mutate(est_type = paste(est_type, " (+Weights)"))
  bind_rows(bivar, covs, covs_wt) %>% mutate(section = section_name)
})

outcome_groups <- list(
  Breastfeeding = c('q19', 'q20', 'q24', 'q25', 'q26', 
                    'q27', 'q28', 'q30', 'q31', 'breastfeed_knowledge'),
  
  Pregnancy = c('q34', 'q36', 'q37', 'q40', 'q42'),
  
  Food_Supplements = c('q43wt', 'q46wt', 'q47', 'q44'),
  
  WASH = c('q65', 'q66', 'wash_index')
  
)

variable_labels <- tibble(
  outcome = unlist(outcome_groups),  
  label = unlist(c(
    c(
      "Children Alive",
      "Has Health Card",
      "Time to First Breastfeed",
      "Non-Breastfeeding in First 2 Days",
      "Still Breastfeeding",
      "Breastfed Times in Last 24 Hrs",
      "Age When Supplemented",
      "Planned Months Exclusive BF",
      "Planned Months Supplementary BF",
      "Breastfeeding Knowledge"
    ),
    
    c(
      "Did you ever receive training on good pregnancy habits in last year?",
      "Months Pregnant at First ANC Visit",
      "Number of ANC Visits During Pregnancy",
      "Took Nutritional Supplements During Pregnancy",
      "Days Took Tablets or Syrup During Pregnancy"
    ),
    
    c(
      "Respondent Food Consumption Score (FCS)",
      "Child Food Consumption Score (FCS)",
      "Nutrition Supplements (Adults)",
      "Nutrition Supplements (Children)"
    ),
    
    c(
      "Wash Hands: Cooking",
      "Wash Hands: Eating",
      "Wash Hands: Index"
    )
  )))

results_all <- results_all %>%
  left_join(variable_labels, by = "outcome") %>%
  mutate(label_wrapped = stringr::str_wrap(label, width = 25))

plot_by_section("Breastfeeding", results_all)

plot_by_section("Pregnancy", results_all)

plot_by_section("Food_Supplements", results_all)

plot_by_section('WASH', results_all)

```

```{r}
# ------------------------------------------------------------------------------
# Helper function --------------------------------------------------------------
# ------------------------------------------------------------------------------

run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- se(model)["treat"]
      tibble::tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = nobs(model)
      )
    } else {
      tibble::tibble(
        outcome = outcome, estimate = NA_real_, std_error = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
        control_mean = NA_real_, N = NA_integer_
      )
    }
  })
}

results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>% 
             dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>% dplyr::mutate(section = section_name)
}) %>%
  dplyr::left_join(variable_labels, by = "outcome") %>%
  dplyr::mutate(label_wrapped = stringr::str_wrap(label, width = 25))

bd_table_bivar    <- format_thp_table(results_all, "^bivariate$")
bd_table_controls <- format_thp_table(results_all, "^controls$")
bd_table_weighted <- format_thp_table(results_all, "controls\\s+\\(\\+weights\\)")

print(bd_table_bivar, n = Inf, width = Inf)
print(bd_table_controls, n = Inf, width = Inf)
print(bd_table_weighted, n = Inf, width = Inf)

write.csv(bd_table_bivar,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_bivar.csv")
write.csv(bd_table_controls,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_controls.csv")
write.csv(bd_table_weighted,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_weighted.csv")
```

```{r bangladesh-revised}

# -----------------------------------------------
# Load & build variables
# -----------------------------------------------
thp_raw <- haven::read_sav('F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_Clean_Data.sav')

thp_clean <- thp_raw %>%
  dplyr::mutate(across(where(is.numeric), ~ if_else(. %in% c(99, 999, 9999), NA_real_, .))) %>%
  dplyr::mutate(
    yes_e3 = if_else(e3 == 1, 1, 0),

    q20 = if_else(q20 %in% c('1', '2'), 1, 0),
    q24 = if_else(q24 == 999, NA_real_, q24),
    q28 = if_else(q28 %in% c(96, 99, 100), NA, q28),

    correct_q32 = if_else(q32 == 3, 1, 0),
    correct_q33 = if_else(q33 == 6, 1, 0),
    breastfeed_knowledge = (correct_q32 + correct_q33) / 2,

    wash_index = (q65 + q66) / 12
  )

# FCS-style weighted sums
thp_clean$q43wt <- 2 * rowSums(thp_clean[, c("q43_03", "q43_04", "q43_01", "q43_02")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q43_05", "q43_06")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q43_07", "q43_08", "q43_09")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q43_012","q43_013","q43_014","q43_015","q43_016",
                              "q43_017","q43_018","q43_019","q43_020","q43_021")], na.rm = TRUE)

thp_clean$q46wt <- 2 * rowSums(thp_clean[, c("q46_03", "q46_04", "q46_01", "q46_02")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q46_05", "q46_06")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q46_07", "q46_08", "q46_09")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q46_012","q46_013","q46_014","q46_015","q46_016",
                              "q46_017","q46_018","q46_019","q46_020","q46_021")], na.rm = TRUE)

# Wealth & other recodes
thp_clean <- thp_clean %>%
  dplyr::mutate(across(c(q2, q12, q13), ~ as.factor(.))) %>%
  dplyr::mutate(q20 = if_else(q20 %in% c(1, 2), 1, 0)) %>%
  dplyr::mutate(
    q5   = if_else(q5 %in% c('96','99'), NA, as.numeric(q5)),
    roof = if_else(q9 > 3, 1, 0),
    floor= if_else(q8 > 4, 1, 0),
    water= if_else(q10 %in% c(6,3,5), 1, 0),
    phone= if_else(q11 %in% c(2,3,4,5), 1, 0)
  ) %>%
  dplyr::mutate(wlth = prcomp(dplyr::select(., roof, floor, water, phone))$x[,1]) %>%
  dplyr::mutate(
    q18_w = suppressWarnings(as.numeric(q18))/4,
    q21_w = suppressWarnings(as.numeric(as.Date(b13) - as.Date(q21)))/7
  ) %>%
  dplyr::mutate(treat = yes_e3)

# -----------------------------------------------
# Specs
# -----------------------------------------------
control_vars <- c('q2','q3','q5','q6','q7','q12','q13','wlth','q16')

outcome_groups <- list(
  Breastfeeding    = c('q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge'),
  Pregnancy        = c('q34','q36','q37','q40','q42'),
  Food_Supplements = c('q43wt','q46wt','q44','q47'),
  WASH             = c('q65','q66','wash_index')
)

# (optional pretty labels; remove this block if you prefer raw names in output)
variable_labels <- tibble::tibble(
  outcome = unlist(outcome_groups),
  label   = unlist(c(
    c("Children Alive","Has Health Card","Time to First Breastfeed","Non-Breastfeeding in First 2 Days",
      "Still Breastfeeding","Breastfed Times in Last 24 Hrs","Age When Supplemented",
      "Planned Months Exclusive BF","Planned Months Supplementary BF","Breastfeeding Knowledge"),
    c("Training on Pregnancy Habits (12m)","Months Pregnant at 1st ANC","Number of ANC Visits",
      "Supplements During Pregnancy","Days Took Tablets/Syrup"),
    c("Respondent FCS","Child FCS","Nutrition Supplements (Adults)","Nutrition Supplements (Children)"),
    c("Wash Hands: Cooking","Wash Hands: Eating","Wash Hands: Index")
  ))
)

# -----------------------------------------------
# Propensity, trimming, stabilized IPW
# -----------------------------------------------
fml <- as.formula(paste0('treat ~ ', paste(control_vars, collapse = " + ")))
pred_treat <- glm(fml, data = thp_clean, family = binomial(link = 'logit'))
thp_clean$prob <- predict(pred_treat, type = 'response')

# Trim to overlap like Mexico
thp_final <- thp_clean %>%
  dplyr::filter(prob > 0.10, prob < 0.75)

# Stabilized IPW: Hájek-style
pbar <- mean(thp_final$treat)
thp_final <- thp_final %>%
  dplyr::mutate(
    weight = if_else(treat == 1, pbar / prob, (1 - pbar) / (1 - prob))
  )

# -----------------------------------------------
# Control-relative standardization (Z on control)
# -----------------------------------------------
vars_to_std <- unique(unlist(outcome_groups))

ctrl_means <- thp_final %>%
  dplyr::filter(treat == 0) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(vars_to_std), ~ mean(., na.rm = TRUE))) %>%
  as.list() |> unlist()

ctrl_sds <- thp_final %>%
  dplyr::filter(treat == 0) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(vars_to_std), ~ sd(., na.rm = TRUE))) %>%
  as.list() |> unlist()
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

thp_scale <- thp_final %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()]
  ))

# -----------------------------------------------
# Estimation helpers (same three specs)
# -----------------------------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(fixest::feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- fixest::se(model)["treat"]
      tibble::tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = stats::nobs(model)
      )
    } else {
      tibble::tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
                     ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
                     control_mean = NA_real_, N = NA_integer_)
    }
  })
}

results_all <- purrr::imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>%
             dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>%
    dplyr::mutate(section = section_name)
})

# (optional) attach labels for presentation
results_all <- results_all %>%
  dplyr::left_join(variable_labels, by = "outcome") %>%
  dplyr::mutate(label_wrapped = stringr::str_wrap(label, width = 25))

# -----------------------------------------------
# Table formatter (control_mean ~ 0 by construction)
# -----------------------------------------------
format_thp_table <- function(results_all, type_pattern){
  df <- results_all %>% dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE))
  if (!"N" %in% names(df)) df <- dplyr::mutate(df, N = NA_integer_)
  df %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 2),
      se        = round(se, 2),
      ci        = paste0("[", round(conf.low, 2), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 2),
      outcome   = dplyr::if_else(!is.na(label), label, outcome)
    ) %>%
    dplyr::select(section, outcome, att, se, ci, n, ctrl_mean)
}

bd_table_bivar    <- format_thp_table(results_all, "^bivariate$")
bd_table_controls <- format_thp_table(results_all, "^controls$")
bd_table_weighted <- format_thp_table(results_all, "controls\\s+\\(\\+weights\\)")

print(bd_table_bivar, n = Inf, width = Inf)
print(bd_table_controls, n = Inf, width = Inf)
print(bd_table_weighted, n = Inf, width = Inf)

# Save
write.csv(bd_table_bivar,    "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_bivar_revised.csv",    row.names = FALSE)
write.csv(bd_table_controls, "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_controls_revised.csv", row.names = FALSE)
write.csv(bd_table_weighted, "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_weighted_revised.csv", row.names = FALSE)
```

```{r bagladesh-revised-all-regressors}
# ------------------------------------------------------------
# 0) Define index components exactly as used to build the indexes
#     (keep only those that exist in the data to avoid errors)
# ------------------------------------------------------------
components_q43 <- c(
  "q43_03","q43_04","q43_01","q43_02",    # weight 2
  "q43_05","q43_06",                      # weight 1
  "q43_07","q43_08","q43_09",             # weight 4
  "q43_012","q43_013","q43_014","q43_015","q43_016",
  "q43_017","q43_018","q43_019","q43_020","q43_021"  # weight 0.5
)
components_q43 <- intersect(components_q43, names(thp_clean))

components_q46 <- c(
  "q46_03","q46_04","q46_01","q46_02",    # weight 2
  "q46_05","q46_06",                      # weight 1
  "q46_07","q46_08","q46_09",             # weight 4
  "q46_012","q46_013","q46_014","q46_015","q46_016",
  "q46_017","q46_018","q46_019","q46_020","q46_021"  # weight 0.5
)
components_q46 <- intersect(components_q46, names(thp_clean))

# Breastfeeding knowledge components
bf_knowledge_components <- c("correct_q32","correct_q33")
bf_knowledge_components <- intersect(bf_knowledge_components, names(thp_clean))

# ------------------------------------------------------------
# 1) Expand outcomes: include indexes AND their components
# ------------------------------------------------------------
outcome_groups <- list(
  Breastfeeding = c(
    'q19','q20','q24','q25','q26','q27','q28','q30','q31',
    'breastfeed_knowledge', bf_knowledge_components
  ),
  Pregnancy = c('q34','q36','q37','q40','q42'),
  Food_Supplements = c('q43wt','q46wt','q44','q47', components_q43, components_q46),
  WASH = c('q65','q66','wash_index')
)

# ------------------------------------------------------------
# 2) Standardize outcomes relative to the CONTROL group (only change)
#     (Keep all other Bangladesh specs as-is)
# ------------------------------------------------------------
vars_to_std <- unique(unlist(outcome_groups))

# If prob not already computed upstream, compute it (keeps your spec otherwise)
if (!("prob" %in% names(thp_clean))) {
  fml <- as.formula(paste0('yes_e3 ~ ', paste(control_vars, collapse = " + ")))
  pred_treat <- glm(fml, data = thp_clean, family = binomial(link = 'logit'))
  thp_clean$prob <- predict(pred_treat, type = 'response')
}

# Build thp_scale from thp_clean (as in your original), but z-score on control stats
ctrl_means <- thp_clean %>%
  dplyr::filter(yes_e3 == 0) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(vars_to_std), ~ mean(., na.rm = TRUE))) %>%
  as.list() |> unlist()

ctrl_sds <- thp_clean %>%
  dplyr::filter(yes_e3 == 0) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(vars_to_std), ~ sd(., na.rm = TRUE))) %>%
  as.list() |> unlist()
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

thp_scale <- thp_clean %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()]
  ))

# Keep your original treatment/weight spec unchanged
thp_scale$treat  <- thp_scale$yes_e3
thp_scale$weight <- ifelse(thp_scale$treat == 1, 1 / thp_scale$prob, 1 / (1 - thp_scale$prob))

# ------------------------------------------------------------
# 3) Estimation helpers (unchanged)
# ------------------------------------------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(fixest::feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- fixest::se(model)["treat"]
      tibble::tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = stats::nobs(model)
      )
    } else {
      tibble::tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
                     ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
                     control_mean = NA_real_, N = NA_integer_)
    }
  })
}

results_all <- purrr::imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>%
             dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>%
    dplyr::mutate(section = section_name)
})

# ------------------------------------------------------------
# 4) Format & export (keeps raw var names; no relabelling)
# ------------------------------------------------------------
format_thp_table <- function(results_all, type_pattern){
  df <- results_all %>% dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE))
  if (!"N" %in% names(df)) df <- dplyr::mutate(df, N = NA_integer_)
  df %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 2),
      se        = round(se, 2),
      ci        = paste0("[", round(conf.low, 2), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 2)
    ) %>%
    dplyr::select(section, outcome, att, se, ci, n, ctrl_mean)
}

bd_table_bivar_all_regressors    <- format_thp_table(results_all, "^bivariate$")
bd_table_controls_all_regressors <- format_thp_table(results_all, "^controls$")
bd_table_weighted_all_regressors <- format_thp_table(results_all, "controls\\s+\\(\\+weights\\)")

print(bd_table_bivar_all_regressors,    n = Inf, width = Inf)
print(bd_table_controls_all_regressors, n = Inf, width = Inf)
print(bd_table_weighted_all_regressors, n = Inf, width = Inf)

write.csv(bd_table_bivar_all_regressors, "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_bivar_all_regressors.csv",    row.names = FALSE)
write.csv(bd_table_controls_all_regressors, "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_controls_all_regressors.csv", row.names = FALSE)
write.csv(bd_table_weighted_all_regressors, "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/bd_weighted_all_regressors.csv", row.names = FALSE)
```

# THP Zambia

```{r}
thp_raw <- read.csv('C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
                    header = TRUE)

thp_clean <- thp_raw %>%
  mutate(across(where(is.numeric), ~ if_else(. %in% c(99, 999, 9999), NA_real_, .))) %>%
  mutate(
    treat = if_else(e3 == 1, 1, 0),
    q20   = if_else(q20 %in% c('1', '2'), 1, 0),
    q24   = if_else(q24 == 999, NA_real_, q24),
    q28   = if_else(q28 %in% c(96, 99, 100), NA, q28),
    correct_Q32 = if_else(q32 == 3, 1, 0),
    correct_Q33 = if_else(q33 == 6, 1, 0),
    breastfeed_knowledge = (correct_Q32 + correct_Q33) / 2,
    norms_index = (q64 + q65 + (4 - q66)) / 12,
    wash_index  = (q68 + q69) / 12
  )

# Weighted FCS (your original mapping)
thp_clean$q45wt <- 2 * rowSums(thp_clean[, c("q453","q454","q451","q452")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q455","q456")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q457","q458","q459")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4512","q4513","q4514","q4515","q4516",
                              "q4517","q4518","q4519","q4520","q4521")], na.rm = TRUE)

# Weighted FCS (your original mapping)
thp_clean$q48wt <- 2 * rowSums(thp_clean[, c("q483","q484","q481","q482")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q485","q486")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q487","q488","q489")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4812","q4813","q4814","q4815","q4816",
                              "q4817", "q4819","q4820","q4821")], na.rm = TRUE)

# Anthropometrics (assumes anthro_zscores() is available)
anthro_data <- thp_clean %>%
  select(case_id, q22, a11, a33, a22, a44)

anthro_data$muac_status <- cut(anthro_data$a44,
                               breaks = c(-Inf, 11.5, 12.5, Inf),
                               labels = c("SAM", "MAM", "Normal"))

results <- anthro_zscores(
  sex = anthro_data$q22,
  age = anthro_data$a11,
  is_age_in_month = TRUE,
  weight = anthro_data$a33,
  lenhei = anthro_data$a22,
  armc = anthro_data$a44,
  measure = "h"
)

results <- results %>%
  select(zwei, zac, zwfl, zlen) %>%
  rename(
    weight_for_age     = zwei,
    arm_for_age        = zac,
    weight_for_height  = zwfl,
    height_for_age     = zlen
  ) %>%
  mutate(
    muac_status = anthro_data$muac_status,
    muac_normal = as.numeric(muac_status == "Normal")
  )

thp_clean <- thp_clean %>% bind_cols(results)

# Keep only variables you listed
histogram_vars <- c(
  'q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge',
  'q34','q36','q37','q40','q42',
  'q45wt', 'q48wt', 'q46','q49',
  'weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal',
  'q64','q65','q66','norms_index',
  'q68','q69','wash_index',
  'tz1','tz3','tz4','tz6',
  'tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25',
  'tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b',
  'treat'
)
thp_clean <- thp_clean %>% select(all_of(histogram_vars))

# --- 1) Second read (robust recodes for wlth etc.) + safe PCA -----------------
safe_pc1 <- function(df, vars){
  M  <- dplyr::select(df, dplyr::all_of(vars))
  ok <- stats::complete.cases(M)
  out <- rep(NA_real_, nrow(M))
  if (sum(ok) >= 2) {
    pc <- stats::prcomp(M[ok, ], center = TRUE, scale. = TRUE)
    out[ok] <- pc$x[,1]
  }
  out
}

thp_covs <- read.csv('C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
                     header = TRUE) %>%
  mutate(
    # numeric & text paths (handles either file flavor without changing your logic)
    q8_num  = suppressWarnings(as.numeric(q8)),
    q9_num  = suppressWarnings(as.numeric(q9)),
    q10_num = suppressWarnings(as.numeric(q10)),
    q11_num = suppressWarnings(as.numeric(q11)),
    q20_num = suppressWarnings(as.numeric(q20)),
    e3_num  = suppressWarnings(as.numeric(e3)),
    q8_chr  = str_trim(as.character(q8)),
    q9_chr  = str_trim(as.character(q9)),
    q10_chr = str_trim(as.character(q10)),
    q11_chr = str_trim(as.character(q11)),
    q20_chr = str_trim(as.character(q20)),
    e3_chr  = str_trim(as.character(e3)),

    # same q20 rule (codes 1/2 OR seen/seen-not)
    q20 = if_else((!is.na(q20_num) & q20_num %in% c(1,2)) |
                  q20_chr %in% c('Yes, seen','Yes, but not seen',' Yes, seen',' Yes, but not seen'), 1L, 0L),

    # asset dummies (codes OR labels, same categories you used)
    roof = case_when(
      (!is.na(q9_num) & q9_num %in% c(9,5,4)) |
      q9_chr %in% c('9 Cement','5 Metal / tin','4 Wood planks',' 9 Cement',' 5 Metal / tin',' 4 Wood planks') ~ 1L,
      (!is.na(q9_num) & q9_num %in% c(2,3)) |
      q9_chr %in% c('2 Thatch / palm leaf / straw','3 Palm tree / bamboo','other',' other',
                    ' 2 Thatch / palm leaf / straw',' 3 Palm tree / bamboo') ~ 0L,
      TRUE ~ NA_integer_
    ),
    floor = case_when(
      (!is.na(q8_num) & q8_num == 7) | q8_chr %in% c('7 Cement',' 7 Cement') ~ 1L,
      is.na(q8) ~ NA_integer_,
      TRUE ~ 0L
    ),
    water = case_when(
      (!is.na(q10_num) & q10_num %in% c(6,3,5)) |
      q10_chr %in% c('6 Protected well','3 Piped to neighbor','5 Tube well or borehole',
                     ' 6 Protected well',' 3 Piped to neighbor',' 5 Tube well or borehole') ~ 1L,
      is.na(q10) ~ NA_integer_,
      TRUE ~ 0L
    ),
    phone = case_when(
      (!is.na(q11_num) & q11_num == 1) | q11_chr %in% c('1 Yes',' 1 Yes') ~ 1L,
      is.na(q11) ~ NA_integer_,
      TRUE ~ 0L
    ),

    q18_w = suppressWarnings(as.numeric(q18))/4,
    q21_w = as.numeric(as.Date(submissiondate) - as.Date(q21))/7,
    q17   = if_else(str_trim(as.character(q17)) == 'Not Applicable',
                    NA_real_, suppressWarnings(as.numeric(q17))),
    yes_e3 = if_else((!is.na(e3_num) & e3_num == 1) | e3_chr %in% c('1','1 Yes',' 1 Yes'), 1L, 0L)
  ) %>%
  mutate(
    wlth = safe_pc1(., c('roof','floor','water','phone'))
  ) %>%
  select(q2, q3, q5, q6, q7, q12, q13, wlth, q16, q17, q18_w, q21_w, yes_e3)

# bind covariates to outcomes/anthro
thp_clean <- cbind.data.frame(thp_covs, thp_clean)
if (!"treat" %in% names(thp_clean)) thp_clean$treat <- thp_clean$yes_e3

# --- 2) Balance, propensity, trimming (unchanged thresholds) ------------------
control_vars <- c('q2','q3','q6','q7','q12','q13','wlth','q16')
fml <- as.formula(paste0('treat ~ ', paste(control_vars, collapse = " + ")))

# (Optional) vtable::st(...) balance tables can remain as you had them

pred_treat <- glm(fml, data = thp_clean, family = binomial(link = 'logit'))
thp_clean$prob <- predict(pred_treat, type = 'response')

thp_clean <- thp_clean %>%
  filter(prob > .3 & prob < .5)

# --- 3) Scale outcomes (fix parentheses!) & stabilized weights ----------------
thp_scale <- thp_clean %>%
  mutate(across(
    c(
      # Breastfeeding
      'q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge',
      # Pregnancy
      'q34','q36','q37','q40','q42',
      # Food & supplements
      'q45wt', 'q48wt', 'q46','q49',
      # Anthropometry
      'weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal',
      # Norms & WASH
      'q64','q65','q66','norms_index','q68','q69','wash_index',
      # Participation / Knowledge / Behaviors
      'tz1','tz3','tz4','tz6','tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25',
      'tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b'
    ),
    ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
  )) %>%
  mutate(treat = yes_e3) %>%
  mutate(
    weight = if_else(treat == 1, mean(treat, na.rm = TRUE) / prob,
                                (1 - mean(treat, na.rm = TRUE)) / (1 - prob))
  )

# --- 4) Outcomes, labels, models, and HKI-style tables ------------------------
outcome_groups <- list(
  Breastfeeding     = c('q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge'),
  Pregnancy         = c('q34','q36','q37','q40','q42'),
  Food_Supplements  = c('q45wt', 'q48wt', 'q46','q49'),
  Anthropometry     = c('weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal'),
  Norms_Stigmas     = c('q64','q65','q66','norms_index'),
  WASH              = c('q68','q69','wash_index'),        # <- fix: "wash_index" (not "index")
  Participation     = c('tz1','tz3','tz4','tz6'),
  Knowledge_Skills  = c('tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25'),
  Behaviors         = c('tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b')
)

variable_labels <- tibble(
  outcome = unlist(outcome_groups),
  label = unlist(c(
    c(
      "Children Alive","Has Health Card","Time to First Breastfeed","Non-Breastfeeding in First 2 Days",
      "Still Breastfeeding","Breastfed Times in Last 24 Hrs","Age When Supplemented",
      "Planned Months Exclusive BF","Planned Months Supplementary BF","Breastfeeding Knowledge"
    ),
    c(
      "Did you ever receive training on good pregnancy habits in last year?",
      "Months Pregnant at First ANC Visit","Number of ANC Visits During Pregnancy",
      "Took Nutritional Supplements During Pregnancy","Days Took Tablets or Syrup During Pregnancy"
    ),
    c("Food Consumption Score (FCS)", "Respondent Food Consumption Score (FCS)", 
      "Nutrition Supplements (Adults)","Nutrition Supplements (Children)"),
    c("Weight for Age","MUAC for Age","Weight for Height","Height for Age","Share with Normal MUAC"),
    c("Agree: Supplements Means Failure as Caregiver","Agree: Supplements Cause Large Baby & Painful Birth",
      "Confidence in Positioning & Attaching Baby for Breastfeeding","Norms & Stigma Index"),
    c("Wash Hands: Cooking","Wash Hands: Eating","Wash Hands: Index"),
    c("Belongs to Mother Support Group","MSG Meets Regularly","Frequency of MSG Meetings","Times Attended MSG"),
    c("Has Wooden Rack for Drying Dishes","Has Access to Latrine","Tippy-Tap Close to Latrine",
      "Village Certified Open Defecation Free","Confidence in Drying Foods","Believes Vaccines Are Good",
      "Confidence in Preserving Supplements","Knows How to Prepare Fortified Porridge","Boils Water for Drinking"),
    c("Drying Fruits, Vegetables, Meat","Turning Soya/Groundnuts into Milk/Biscuits",
      "Prepared Nutritious Food at Home","None of the Above","Received Packets of Seeds","Has Home/Keyhole Garden",
      "Confidence Growing Vegetables","Garden Meets Family Nutrition Needs")
  ))
)

run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- se(model)["treat"]
      tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = nobs(model)
      )
    } else {
      tibble(
        outcome = outcome, estimate = NA_real_, std_error = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
        control_mean = NA_real_, N = NA_integer_
      )
    }
  })
}

results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>%
             mutate(est_type = paste(est_type, " (+Weights)"))
  bind_rows(bivar, covs, covs_wt) %>% mutate(section = section_name)
}) %>%
  left_join(variable_labels, by = "outcome") %>%
  mutate(label_wrapped = stringr::str_wrap(label, width = 25))

# HKI-style formatter
format_hki_table <- function(results_all, type_pattern){
  results_all %>%
    dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE)) %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 3),
      se        = round(se, 3),
      ci        = paste0("[", round(conf.low, 3), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 3),
      outcome   = dplyr::if_else(!is.na(label), label, outcome)
    ) %>%
    dplyr::select(outcome, att, se, ci, n, ctrl_mean)
}

zmb_table_bivar    <- format_hki_table(results_all, "^bivariate$")
zmb_table_controls <- format_hki_table(results_all, "^controls$")
zmb_table_weighted <- format_hki_table(results_all, "controls\\s+\\(\\+weights\\)")

# Inspect in console
print(zmb_table_bivar, n = Inf, width = Inf)
print(zmb_table_controls, n = Inf, width = Inf)
print(zmb_table_weighted, n = 50, width = Inf)

write.csv(zmb_table_bivar,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_bivar.csv")
write.csv(zmb_table_controls,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_controls.csv")
write.csv(zmb_table_weighted,
          "C:/Users/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_weighted.csv")
```

```{r zambia-revised}
# -----------------------------
# 1) Load raw & primary recodes
# -----------------------------
thp_raw <- read.csv(
  'F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
  header = TRUE
)

thp_clean <- thp_raw %>%
  mutate(across(where(is.numeric), ~ if_else(. %in% c(99, 999, 9999), NA_real_, .))) %>%
  mutate(
    treat  = if_else(e3 == 1, 1, 0),
    q20    = if_else(q20 %in% c('1','2'), 1, 0),
    q24    = if_else(q24 == 999, NA_real_, q24),
    q28    = if_else(q28 %in% c(96, 99, 100), NA, q28),
    correct_Q32 = if_else(q32 == 3, 1, 0),
    correct_Q33 = if_else(q33 == 6, 1, 0),
    breastfeed_knowledge = (correct_Q32 + correct_Q33) / 2,
    norms_index = (q64 + q65 + (4 - q66)) / 12,
    wash_index  = (q68 + q69) / 12
  )

# -----------------------------
# 2) FCS-style composite scores
# -----------------------------
thp_clean$q45wt <- 2 * rowSums(thp_clean[, c("q453","q454","q451","q452")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q455","q456")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q457","q458","q459")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4512","q4513","q4514","q4515","q4516",
                              "q4517","q4518","q4519","q4520","q4521")], na.rm = TRUE)

thp_clean$q48wt <- 2 * rowSums(thp_clean[, c("q483","q484","q481","q482")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q485","q486")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q487","q488","q489")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4812","q4813","q4814","q4815","q4816",
                              "q4817","q4819","q4820","q4821")], na.rm = TRUE)

# -----------------------------
# 3) Anthropometrics (assumes anthro_zscores() exists in your env)
# -----------------------------
anthro_data <- thp_clean %>%
  select(case_id, q22, a11, a33, a22, a44)

anthro_data$muac_status <- cut(anthro_data$a44,
                               breaks = c(-Inf, 11.5, 12.5, Inf),
                               labels = c("SAM", "MAM", "Normal"))

results <- anthro_zscores(
  sex = anthro_data$q22,
  age = anthro_data$a11,
  is_age_in_month = TRUE,
  weight = anthro_data$a33,
  lenhei = anthro_data$a22,
  armc = anthro_data$a44,
  measure = "h"
)

results <- results %>%
  select(zwei, zac, zwfl, zlen) %>%
  rename(
    weight_for_age     = zwei,
    arm_for_age        = zac,
    weight_for_height  = zwfl,
    height_for_age     = zlen
  ) %>%
  mutate(
    muac_status = anthro_data$muac_status,
    muac_normal = as.numeric(muac_status == "Normal")
  )

thp_clean <- thp_clean %>% bind_cols(results)

# Keep the variables of interest (as in your list)
histogram_vars <- c(
  'q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge',
  'q34','q36','q37','q40','q42',
  'q45wt','q48wt','q46','q49',
  'weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal',
  'q64','q65','q66','norms_index',
  'q68','q69','wash_index',
  'tz1','tz3','tz4','tz6',
  'tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25',
  'tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b',
  'treat'
)
thp_clean <- thp_clean %>% select(all_of(histogram_vars))

# ------------------------------------------------------
# 4) Build covariates & wealth index from the original file
# ------------------------------------------------------
safe_pc1 <- function(df, vars){
  M  <- dplyr::select(df, dplyr::all_of(vars))
  ok <- stats::complete.cases(M)
  out <- rep(NA_real_, nrow(M))
  if (sum(ok) >= 2) {
    pc <- stats::prcomp(M[ok, ], center = TRUE, scale. = TRUE)
    out[ok] <- pc$x[,1]
  }
  out
}

thp_covs <- read.csv(
  'F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
  header = TRUE
) %>%
  mutate(
    q8_num  = suppressWarnings(as.numeric(q8)),
    q9_num  = suppressWarnings(as.numeric(q9)),
    q10_num = suppressWarnings(as.numeric(q10)),
    q11_num = suppressWarnings(as.numeric(q11)),
    q20_num = suppressWarnings(as.numeric(q20)),
    e3_num  = suppressWarnings(as.numeric(e3)),

    q8_chr  = str_trim(as.character(q8)),
    q9_chr  = str_trim(as.character(q9)),
    q10_chr = str_trim(as.character(q10)),
    q11_chr = str_trim(as.character(q11)),
    q20_chr = str_trim(as.character(q20)),
    e3_chr  = str_trim(as.character(e3)),

    q20 = if_else((!is.na(q20_num) & q20_num %in% c(1,2)) |
                  q20_chr %in% c('Yes, seen','Yes, but not seen',' Yes, seen',' Yes, but not seen'), 1L, 0L),

    roof = case_when(
      (!is.na(q9_num) & q9_num %in% c(9,5,4)) |
      q9_chr %in% c('9 Cement','5 Metal / tin','4 Wood planks',' 9 Cement',' 5 Metal / tin',' 4 Wood planks') ~ 1L,
      (!is.na(q9_num) & q9_num %in% c(2,3)) |
      q9_chr %in% c('2 Thatch / palm leaf / straw','3 Palm tree / bamboo','other',' other',
                    ' 2 Thatch / palm leaf / straw',' 3 Palm tree / bamboo') ~ 0L,
      TRUE ~ NA_integer_
    ),
    floor = case_when(
      (!is.na(q8_num) & q8_num == 7) | q8_chr %in% c('7 Cement',' 7 Cement') ~ 1L,
      is.na(q8) ~ NA_integer_,
      TRUE ~ 0L
    ),
    water = case_when(
      (!is.na(q10_num) & q10_num %in% c(6,3,5)) |
      q10_chr %in% c('6 Protected well','3 Piped to neighbor','5 Tube well or borehole',
                     ' 6 Protected well',' 3 Piped to neighbor',' 5 Tube well or borehole') ~ 1L,
      is.na(q10) ~ NA_integer_,
      TRUE ~ 0L
    ),
    phone = case_when(
      (!is.na(q11_num) & q11_num == 1) | q11_chr %in% c('1 Yes',' 1 Yes') ~ 1L,
      is.na(q11) ~ NA_integer_,
      TRUE ~ 0L
    ),

    q18_w = suppressWarnings(as.numeric(q18))/4,
    q21_w = as.numeric(as.Date(submissiondate) - as.Date(q21))/7,
    q17   = if_else(str_trim(as.character(q17)) == 'Not Applicable',
                    NA_real_, suppressWarnings(as.numeric(q17))),
    yes_e3 = if_else((!is.na(e3_num) & e3_num == 1) | e3_chr %in% c('1','1 Yes',' 1 Yes'), 1L, 0L)
  ) %>%
  mutate(
    wlth = safe_pc1(., c('roof','floor','water','phone'))
  ) %>%
  select(q2, q3, q5, q6, q7, q12, q13, wlth, q16, q17, q18_w, q21_w, yes_e3)

# Attach covariates to outcomes/anthro
thp_clean <- cbind.data.frame(thp_covs, thp_clean)
if (!"treat" %in% names(thp_clean)) thp_clean$treat <- thp_clean$yes_e3

# -----------------------------
# 5) Propensity & trimming
# -----------------------------
control_vars <- c('q2','q3','q6','q7','q12','q13','wlth','q16')
fml <- as.formula(paste0('treat ~ ', paste(control_vars, collapse = " + ")))

pred_treat <- glm(fml, data = thp_clean, family = binomial(link = 'logit'))
thp_clean$prob <- predict(pred_treat, type = 'response')

# Your original overlap window
thp_clean <- thp_clean %>% filter(prob > .3 & prob < .5)

# ----------------------------------------------------
# 6) Control-relative standardization (MAIN REVISION)
# ----------------------------------------------------
outcome_groups <- list(
  Breastfeeding     = c('q19','q20','q24','q25','q26','q27','q28','q30','q31','breastfeed_knowledge'),
  Pregnancy         = c('q34','q36','q37','q40','q42'),
  Food_Supplements  = c('q45wt','q48wt','q46','q49'),
  Anthropometry     = c('weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal'),
  Norms_Stigmas     = c('q64','q65','q66','norms_index'),
  WASH              = c('q68','q69','wash_index'),
  Participation     = c('tz1','tz3','tz4','tz6'),
  Knowledge_Skills  = c('tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25'),
  Behaviors         = c('tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b')
)

vars_to_std <- unique(unlist(outcome_groups))

ctrl_means <- thp_clean %>%
  filter(treat == 0) %>%
  summarise(across(all_of(vars_to_std), ~ mean(., na.rm = TRUE))) %>%
  as.list() |> unlist()

ctrl_sds <- thp_clean %>%
  filter(treat == 0) %>%
  summarise(across(all_of(vars_to_std), ~ sd(., na.rm = TRUE))) %>%
  as.list() |> unlist()
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

thp_scale <- thp_clean %>%
  mutate(across(
    all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()]
  )) %>%
  mutate(
    treat  = treat,
    weight = if_else(
      treat == 1,
      mean(treat, na.rm = TRUE) / prob,
      (1 - mean(treat, na.rm = TRUE)) / (1 - prob)
    )
  )

# -----------------------------
# 7) Labels (optional, for display)
# -----------------------------
variable_labels <- tibble(
  outcome = unlist(outcome_groups),
  label = unlist(c(
    c(
      "Children Alive","Has Health Card","Time to First Breastfeed","Non-Breastfeeding in First 2 Days",
      "Still Breastfeeding","Breastfed Times in Last 24 Hrs","Age When Supplemented",
      "Planned Months Exclusive BF","Planned Months Supplementary BF","Breastfeeding Knowledge"
    ),
    c(
      "Did you ever receive training on good pregnancy habits in last year?",
      "Months Pregnant at First ANC Visit","Number of ANC Visits During Pregnancy",
      "Took Nutritional Supplements During Pregnancy","Days Took Tablets or Syrup During Pregnancy"
    ),
    c("Food Consumption Score (FCS)", "Respondent Food Consumption Score (FCS)",
      "Nutrition Supplements (Adults)","Nutrition Supplements (Children)"),
    c("Weight for Age","MUAC for Age","Weight for Height","Height for Age","Share with Normal MUAC"),
    c("Agree: Supplements Means Failure as Caregiver","Agree: Supplements Cause Large Baby & Painful Birth",
      "Confidence in Positioning & Attaching Baby for Breastfeeding","Norms & Stigma Index"),
    c("Wash Hands: Cooking","Wash Hands: Eating","Wash Hands: Index"),
    c("Belongs to Mother Support Group","MSG Meets Regularly","Frequency of MSG Meetings","Times Attended MSG"),
    c("Has Wooden Rack for Drying Dishes","Has Access to Latrine","Tippy-Tap Close to Latrine",
      "Village Certified Open Defecation Free","Confidence in Drying Foods","Believes Vaccines Are Good",
      "Confidence in Preserving Supplements","Knows How to Prepare Fortified Porridge","Boils Water for Drinking"),
    c("Drying Fruits, Vegetables, Meat","Turning Soya/Groundnuts into Milk/Biscuits",
      "Prepared Nutritious Food at Home","None of the Above","Received Packets of Seeds","Has Home/Keyhole Garden",
      "Confidence Growing Vegetables","Garden Meets Family Nutrition Needs")
  ))
)

# -----------------------------
# 8) Estimation & tables
# -----------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- se(model)["treat"]
      tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = nobs(model)
      )
    } else {
      tibble(
        outcome = outcome, estimate = NA_real_, std_error = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
        control_mean = NA_real_, N = NA_integer_
      )
    }
  })
}

control_vars <- c('q2','q3','q6','q7','q12','q13','wlth','q16')

results_all <- imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>% mutate(est_type = paste(est_type, " (+Weights)"))
  bind_rows(bivar, covs, covs_wt) %>% mutate(section = section_name)
}) %>%
  left_join(variable_labels, by = "outcome") %>%
  mutate(label_wrapped = stringr::str_wrap(label, width = 25))

format_hki_table <- function(results_all, type_pattern){
  results_all %>%
    dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE)) %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 3),
      se        = round(se, 3),
      ci        = paste0("[", round(conf.low, 3), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 3),
      outcome   = dplyr::if_else(!is.na(label), label, outcome)
    ) %>%
    dplyr::select(outcome, att, se, ci, n, ctrl_mean)
}

zmb_table_bivar    <- format_hki_table(results_all, "^bivariate$")
zmb_table_controls <- format_hki_table(results_all, "^controls$")
zmb_table_weighted <- format_hki_table(results_all, "controls\\s+\\(\\+weights\\)")

print(zmb_table_bivar, n = Inf, width = Inf)
print(zmb_table_controls, n = Inf, width = Inf)
print(zmb_table_weighted, n = 50, width = Inf)

# -----------------------------
# 9) Save CSVs
# -----------------------------
write.csv(zmb_table_bivar,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_bivar_revised.csv",
          row.names = FALSE)
write.csv(zmb_table_controls,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_controls_revised.csv",
          row.names = FALSE)
write.csv(zmb_table_weighted,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_weighted_revised.csv",
          row.names = FALSE)
```

```{r zambia-revised-all-regressors}
# -----------------------------
# 1) Load & primary recodes
# -----------------------------
thp_raw <- read.csv(
  'F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
  header = TRUE
)

thp_clean <- thp_raw %>%
  mutate(across(where(is.numeric), ~ if_else(. %in% c(99, 999, 9999), NA_real_, .))) %>%
  mutate(
    treat  = if_else(e3 == 1, 1, 0),
    q20    = if_else(q20 %in% c('1','2'), 1, 0),
    q24    = if_else(q24 == 999, NA_real_, q24),
    q28    = if_else(q28 %in% c(96, 99, 100), NA, q28),
    correct_Q32 = if_else(q32 == 3, 1, 0),
    correct_Q33 = if_else(q33 == 6, 1, 0),
    breastfeed_knowledge = (correct_Q32 + correct_Q33) / 2,
    norms_index = (q64 + q65 + (4 - q66)) / 12,
    wash_index  = (q68 + q69) / 12
  )

# -----------------------------
# 2) FCS composites
# -----------------------------
thp_clean$q45wt <- 2 * rowSums(thp_clean[, c("q453","q454","q451","q452")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q455","q456")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q457","q458","q459")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4512","q4513","q4514","q4515","q4516",
                              "q4517","q4518","q4519","q4520","q4521")], na.rm = TRUE)

thp_clean$q48wt <- 2 * rowSums(thp_clean[, c("q483","q484","q481","q482")], na.rm = TRUE) +
  1 * rowSums(thp_clean[, c("q485","q486")], na.rm = TRUE) +
  4 * rowSums(thp_clean[, c("q487","q488","q489")], na.rm = TRUE) +
  0.5 * rowSums(thp_clean[, c("q4812","q4813","q4814","q4815","q4816",
                              "q4817","q4819","q4820","q4821")], na.rm = TRUE)

# -----------------------------
# 3) Anthropometrics (assumes anthro_zscores() in env)
# -----------------------------
anthro_data <- thp_clean %>% select(case_id, q22, a11, a33, a22, a44)
anthro_data$muac_status <- cut(anthro_data$a44,
                               breaks = c(-Inf, 11.5, 12.5, Inf),
                               labels = c("SAM", "MAM", "Normal"))

results <- anthro_zscores(
  sex = anthro_data$q22,
  age = anthro_data$a11,
  is_age_in_month = TRUE,
  weight = anthro_data$a33,
  lenhei = anthro_data$a22,
  armc = anthro_data$a44,
  measure = "h"
)

results <- results %>%
  select(zwei, zac, zwfl, zlen) %>%
  rename(
    weight_for_age     = zwei,
    arm_for_age        = zac,
    weight_for_height  = zwfl,
    height_for_age     = zlen
  ) %>%
  mutate(
    muac_status = anthro_data$muac_status,
    muac_normal = as.numeric(muac_status == "Normal")
  )

thp_clean <- thp_clean %>% bind_cols(results)

# ------------------------------------------------------
# 4) Covariates & wealth from original file
# ------------------------------------------------------
safe_pc1 <- function(df, vars){
  M  <- dplyr::select(df, dplyr::all_of(vars))
  ok <- stats::complete.cases(M)
  out <- rep(NA_real_, nrow(M))
  if (sum(ok) >= 2) {
    pc <- stats::prcomp(M[ok, ], center = TRUE, scale. = TRUE)
    out[ok] <- pc$x[,1]
  }
  out
}

thp_covs <- read.csv(
  'F:/cedric.antunes/Documents/Evaluasi/THP_data_input/THP_ZAMBIA_SURVEY_DATA_no_labels.csv',
  header = TRUE
) %>%
  mutate(
    q8_num  = suppressWarnings(as.numeric(q8)),
    q9_num  = suppressWarnings(as.numeric(q9)),
    q10_num = suppressWarnings(as.numeric(q10)),
    q11_num = suppressWarnings(as.numeric(q11)),
    q20_num = suppressWarnings(as.numeric(q20)),
    e3_num  = suppressWarnings(as.numeric(e3)),

    q8_chr  = str_trim(as.character(q8)),
    q9_chr  = str_trim(as.character(q9)),
    q10_chr = str_trim(as.character(q10)),
    q11_chr = str_trim(as.character(q11)),
    q20_chr = str_trim(as.character(q20)),
    e3_chr  = str_trim(as.character(e3)),

    q20 = if_else((!is.na(q20_num) & q20_num %in% c(1,2)) |
                  q20_chr %in% c('Yes, seen','Yes, but not seen',' Yes, seen',' Yes, but not seen'), 1L, 0L),

    roof = case_when(
      (!is.na(q9_num) & q9_num %in% c(9,5,4)) |
      q9_chr %in% c('9 Cement','5 Metal / tin','4 Wood planks',' 9 Cement',' 5 Metal / tin',' 4 Wood planks') ~ 1L,
      (!is.na(q9_num) & q9_num %in% c(2,3)) |
      q9_chr %in% c('2 Thatch / palm leaf / straw','3 Palm tree / bamboo','other',' other',
                    ' 2 Thatch / palm leaf / straw',' 3 Palm tree / bamboo') ~ 0L,
      TRUE ~ NA_integer_
    ),
    floor = case_when(
      (!is.na(q8_num) & q8_num == 7) | q8_chr %in% c('7 Cement',' 7 Cement') ~ 1L,
      is.na(q8) ~ NA_integer_,
      TRUE ~ 0L
    ),
    water = case_when(
      (!is.na(q10_num) & q10_num %in% c(6,3,5)) |
      q10_chr %in% c('6 Protected well','3 Piped to neighbor','5 Tube well or borehole',
                     ' 6 Protected well',' 3 Piped to neighbor',' 5 Tube well or borehole') ~ 1L,
      is.na(q10) ~ NA_integer_,
      TRUE ~ 0L
    ),
    phone = case_when(
      (!is.na(q11_num) & q11_num == 1) | q11_chr %in% c('1 Yes',' 1 Yes') ~ 1L,
      is.na(q11) ~ NA_integer_,
      TRUE ~ 0L
    ),

    q18_w = suppressWarnings(as.numeric(q18))/4,
    q21_w = as.numeric(as.Date(submissiondate) - as.Date(q21))/7,
    q17   = if_else(str_trim(as.character(q17)) == 'Not Applicable',
                    NA_real_, suppressWarnings(as.numeric(q17))),
    yes_e3 = if_else((!is.na(e3_num) & e3_num == 1) | e3_chr %in% c('1','1 Yes',' 1 Yes'), 1L, 0L)
  ) %>%
  mutate(wlth = safe_pc1(., c('roof','floor','water','phone'))) %>%
  select(q2, q3, q5, q6, q7, q12, q13, wlth, q16, q17, q18_w, q21_w, yes_e3)

# Replace overlapping columns with the recoded covariates (keeps thp_covs versions)
overlap <- intersect(names(thp_clean), names(thp_covs))

thp_clean <- dplyr::bind_cols(
  thp_covs,
  thp_clean %>% dplyr::select(-dplyr::all_of(overlap))
)

# ------------------------------------------------------
# 5) INCLUDE all composite components in the dataset
# ------------------------------------------------------
components_q45 <- c("q453","q454","q451","q452","q455","q456","q457","q458","q459",
                    "q4512","q4513","q4514","q4515","q4516","q4517","q4518","q4519","q4520","q4521")
components_q45 <- intersect(components_q45, names(thp_clean))

components_q48 <- c("q483","q484","q481","q482","q485","q486","q487","q488","q489",
                    "q4812","q4813","q4814","q4815","q4816","q4817","q4819","q4820","q4821")
components_q48 <- intersect(components_q48, names(thp_clean))

bf_components <- c("correct_Q32","correct_Q33")
bf_components <- intersect(bf_components, names(thp_clean))

# (No need to reselect columns; the components are already present in thp_clean.)

# -----------------------------
# 6) Propensity & trimming
# -----------------------------
control_vars <- c('q2','q3','q6','q7','q12','q13','wlth','q16')
fml <- as.formula(paste0('treat ~ ', paste(control_vars, collapse = " + ")))

pred_treat <- glm(fml, data = thp_clean, family = binomial(link = 'logit'))
thp_clean$prob <- predict(pred_treat, type = 'response')

# Your overlap window
thp_clean <- thp_clean %>% filter(prob > .3 & prob < .5)

# ----------------------------------------------------
# 7) Outcomes (include indexes + all components)
# ----------------------------------------------------
outcome_groups <- list(
  Breastfeeding     = c('q19','q20','q24','q25','q26','q27','q28','q30','q31',
                        'breastfeed_knowledge', bf_components),
  Pregnancy         = c('q34','q36','q37','q40','q42'),
  Food_Supplements  = c('q45wt','q48wt','q46','q49', components_q45, components_q48),
  Anthropometry     = c('weight_for_age','arm_for_age','weight_for_height','height_for_age','muac_normal'),
  Norms_Stigmas     = c('q64','q65','q66','norms_index'),
  WASH              = c('q68','q69','wash_index'),
  Participation     = c('tz1','tz3','tz4','tz6'),
  Knowledge_Skills  = c('tz17','tz18','tz19','tz20','tz21','tz22','tz23','tz24','tz25'),
  Behaviors         = c('tz140','tz141','tz142','tz143','tz12','tz15','tz16','tz16b')
)

# ----------------------------------------------------
# 8) Control-relative standardization (trimmed sample)
# ----------------------------------------------------
vars_to_std <- unique(unlist(outcome_groups))

ctrl_means <- thp_clean %>%
  filter(treat == 0) %>%
  summarise(across(all_of(vars_to_std), ~ mean(., na.rm = TRUE))) %>%
  as.list() |> unlist()

ctrl_sds <- thp_clean %>%
  filter(treat == 0) %>%
  summarise(across(all_of(vars_to_std), ~ sd(., na.rm = TRUE))) %>%
  as.list() |> unlist()
ctrl_sds[is.na(ctrl_sds) | ctrl_sds == 0] <- 1

thp_scale <- thp_clean %>%
  mutate(across(
    all_of(vars_to_std),
    ~ (. - ctrl_means[cur_column()]) / ctrl_sds[cur_column()]
  )) %>%
  mutate(
    treat  = treat,
    weight = if_else(
      treat == 1,
      mean(treat, na.rm = TRUE) / prob,
      (1 - mean(treat, na.rm = TRUE)) / (1 - prob)
    )
  )

# -----------------------------
# 9) Estimation & tables
# -----------------------------
run_models <- function(outcomes, data, controls = NULL, weights = NULL, model_type = "Bivariate") {
  purrr::map_dfr(outcomes, function(outcome) {
    fml <- as.formula(if (is.null(controls)) paste0(outcome, " ~ treat")
                      else paste0(outcome, " ~ treat + ", paste(controls, collapse = " + ")))
    model <- tryCatch(fixest::feols(fml, data = data, se = "hetero", weights = weights), error = function(e) NULL)

    if (!is.null(model) && "treat" %in% names(coef(model))) {
      est <- coef(model)["treat"]; se_val <- fixest::se(model)["treat"]
      tibble::tibble(
        outcome      = outcome,
        estimate     = est,
        std_error    = se_val,
        ci_low       = est - 1.96 * se_val,
        ci_high      = est + 1.96 * se_val,
        est_type     = model_type,
        control_mean = mean(data[[outcome]][data$treat == 0], na.rm = TRUE),
        N            = stats::nobs(model)
      )
    } else {
      tibble::tibble(outcome = outcome, estimate = NA_real_, std_error = NA_real_,
                     ci_low = NA_real_, ci_high = NA_real_, est_type = model_type,
                     control_mean = NA_real_, N = NA_integer_)
    }
  })
}

results_all <- purrr::imap_dfr(outcome_groups, function(outcomes, section_name) {
  bivar   <- run_models(outcomes, data = thp_scale, model_type = "Bivariate")
  covs    <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls")
  covs_wt <- run_models(outcomes, data = thp_scale, controls = control_vars, model_type = "Controls",
                        weights = thp_scale$weight) %>%
             dplyr::mutate(est_type = paste(est_type, " (+Weights)"))
  dplyr::bind_rows(bivar, covs, covs_wt) %>%
    dplyr::mutate(section = section_name)
})

format_tbl <- function(results_all, type_pattern){
  results_all %>%
    dplyr::filter(grepl(type_pattern, est_type, ignore.case = TRUE)) %>%
    dplyr::rename(conf.low = ci_low, conf.high = ci_high, n = N) %>%
    dplyr::mutate(
      se        = (conf.high - estimate) / 1.96,
      att       = round(estimate, 3),
      se        = round(se, 3),
      ci        = paste0("[", round(conf.low, 3), ", ", round(conf.high, 3), "]"),
      ctrl_mean = round(control_mean, 3)
    ) %>%
    dplyr::select(section, outcome, att, se, ci, n, ctrl_mean)
}

zmb_table_bivar_all_regressors    <- format_tbl(results_all, "^bivariate$")
zmb_table_controls_all_regressors <- format_tbl(results_all, "^controls$")
zmb_table_weighted_all_regressors <- format_tbl(results_all, "controls\\s+\\(\\+weights\\)")

print(zmb_table_bivar_all_regressors, n = Inf, width = Inf)
print(zmb_table_controls_all_regressors, n = Inf, width = Inf)
print(zmb_table_weighted_all_regressors, n = 80, width = Inf)

write.csv(zmb_table_bivar_all_regressors,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_bivar_all_regressors.csv",
          row.names = FALSE)
write.csv(zmb_table_controls_all_regressors,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_controls_all_regressors.csv",
          row.names = FALSE)
write.csv(zmb_table_weighted_all_regressors,
          "F:/cedric.antunes/Documents/Evaluasi/THP_data_output/zamb_weighted_all_regressors.csv",
          row.names = FALSE)
```

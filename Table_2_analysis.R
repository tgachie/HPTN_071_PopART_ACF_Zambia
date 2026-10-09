#===============================================================================#
# Script: 05_table2_analysis.R
#
# Purpose:
# Generate Table 2:
# Factors associated with newly diagnosed TB disease in Round 3
#
# Input:
#   data/acf_regression_r3.csv
#
# Outputs:
#   outputs/tables/Table_2_Associations.csv
#   outputs/tables/Table_2_Associations.docx
#
# Notes:
#   The original primary statistical analysis was conducted using Stata.
#   This R script replicates the Stata survey/cluster-robust logistic regression 
#   models (using finite-sample degrees of freedom corrections) from aggregated data.
#
#   "Unadjusted analysis" reproduces the original Stata analysis in which
#   community was retained in each bivariable model.
#
#   Adjusted analysis includes:
#   community + age + sex + HIV/ART status + household TB contact +
#   residence/mobility status.
#===============================================================================#

#-------------------------------------------------------------------------------#
# 1. Setup & Libraries
#-------------------------------------------------------------------------------#
library(dplyr)
library(readr)
library(tidyr)
library(sandwich)
library(flextable)
library(officer)

# Ensure output directory exists
if (!dir.exists("outputs/tables")) {
  dir.create("outputs/tables", recursive = TRUE)
}

#-------------------------------------------------------------------------------#
# 2. Load Aggregated Round 3 Regression Data
#-------------------------------------------------------------------------------#
dat <- read_csv("data/acf_regression_r3.csv",  show_col_types = FALSE)

#-------------------------------------------------------------------------------#
# 3. Variable Levels and Reference Categories
#-------------------------------------------------------------------------------#
community_numeric_levels <- c(8, 1, 2, 5, 6, 9, 10, 11)
community_labels <- c("Z8", "Z1", "Z2", "Z5", "Z6", "Z9", "Z10", "Z11")
age_levels <- c("15 - 29", "30 - 39", "40 - 49", "50+")
sex_levels <- c("Male", "Female")
hiv_levels <- c("Tested negative by CHiP", "Self-reported living with HIV & NOT on ART",
                "Self-reported living with HIV & on ART", "Tested positive by CHiP",
                "Not tested/indeterminate")
tb_contact_levels <- c("No", "Yes")
movement_levels <- c("Previously resident (same zone)",
                     "Moved within community",
                     "Moved from outside community")

dat <- dat %>%
  mutate(
    com = factor(com,levels = community_numeric_levels, labels = community_labels),
    agecat_2   = factor(agecat_2, levels = age_levels),
    sex        = factor(sex, levels = sex_levels),
    hiv        = factor(hiv, levels = hiv_levels),
    tb_contact = factor(tb_contact, levels = tb_contact_levels),
    movement   = factor(movement, levels = movement_levels),
    
    cluster_id = interaction(com, zone, drop = TRUE)
  )

class(dat$com)
unique(dat$com)

#-------------------------------------------------------------------------------#
# 4. Check Round 3 Totals
#-------------------------------------------------------------------------------#
dat_summary <- dat %>%
  summarise(
    screened = sum(n),
    tb_cases = sum(n_tb),
    yield = 100000 * tb_cases / screened
  )
dat_summary

#===============================================================================#
# SECTION 2: Helper Functions (Stata Replication)
#===============================================================================#

#-------------------------------------------------------------------------------#
# 5. Cluster-Robust Variance Function
#-------------------------------------------------------------------------------#
#
# Stata's cluster-robust variance includes a finite-sample correction.
# Since these data have been aggregated, the correction uses the
# ORIGINAL number of individuals rather than the number of aggregate rows.
#
#-------------------------------------------------------------------------------#

cluster_vcov <- function(model, data) {
  
  # Raw cluster-robust covariance
  V0 <- vcovCL(model, cluster = data$cluster_id, type = "HC0", cadjust = FALSE)
  
  # Number of independent individuals
  N <- sum(data$n)
  
  # Number of clusters
  G <- n_distinct(data$cluster_id)
  
  # Number of fitted parameters
  K <- length(coef(model))
  
  # finite-sample correction
  correction <- (G / (G - 1)) * ((N - 1) / (N - K))
  V0 * correction
}


#-------------------------------------------------------------------------------#
# 6. Overall Wald Test
#-------------------------------------------------------------------------------#
wald_test_terms <- function(model, V,  terms) {
  b <- coef(model)[terms]
  V_sub <- V[terms, terms, drop = FALSE]
  
  W <- as.numeric(t(b) %*% solve(V_sub) %*% b)
  df <- length(b)
  
  pchisq(W, df = df, lower.tail = FALSE)
}


#-------------------------------------------------------------------------------#
# 7. Extract ORs and 95% CIs
#-------------------------------------------------------------------------------#
extract_or <- function(model, V, terms, levels_nonref) {
  
  beta <- coef(model)[terms]
  se <- sqrt(diag(V)[terms])
  
  tibble(level = levels_nonref, 
         estimate = exp(beta), 
         lower = exp(beta - 1.96 * se),
         upper = exp(beta + 1.96 * se))
}


#-------------------------------------------------------------------------------#
# 8. Descriptive Counts and Yield
#-------------------------------------------------------------------------------#
get_descriptive <- function(data, variable,characteristic) {
  
  data %>%
    group_by(level = .data[[variable]]) %>%
    summarise(n_tb = sum(n_tb, na.rm = TRUE),
              N = sum(n, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(yield = 100000 * n_tb / N,
           characteristic = characteristic,
           level = as.character(level)) %>%
    select(characteristic, level, n_tb, N, yield)
}


#===============================================================================#
# SECTION 3: "Unadjusted" Models
#===============================================================================#

# 9. Community
fit_u_com <- glm(cbind(n_tb, n - n_tb) ~ com, 
                 family = binomial(), data = dat)
V_u_com   <- cluster_vcov(fit_u_com, dat)
terms_u_com <- grep("^com", names(coef(fit_u_com)), value = TRUE)
or_u_com  <- extract_or(fit_u_com, V_u_com, terms_u_com, community_levels[-1])
p_u_com   <- wald_test_terms(fit_u_com, V_u_com, terms_u_com)

# 10. Age
fit_u_age <- glm(cbind(n_tb, n - n_tb) ~ com + agecat_2, 
                 family = binomial(), data = dat)
V_u_age   <- cluster_vcov(fit_u_age, dat)
terms_u_age <- grep("^agecat_2", names(coef(fit_u_age)), value = TRUE)
or_u_age  <- extract_or(fit_u_age, V_u_age, terms_u_age, age_levels[-1])
p_u_age   <- wald_test_terms(fit_u_age, V_u_age, terms_u_age)

# 11. Sex
fit_u_sex <- glm(cbind(n_tb, n - n_tb) ~ com + sex, 
                 family = binomial(), data = dat)
V_u_sex   <- cluster_vcov(fit_u_sex, dat)
terms_u_sex <- grep("^sex", names(coef(fit_u_sex)), value = TRUE)
or_u_sex  <- extract_or(fit_u_sex, V_u_sex, terms_u_sex, sex_levels[-1])
p_u_sex   <- wald_test_terms(fit_u_sex, V_u_sex, terms_u_sex)

# 12. HIV/ART Status
fit_u_hiv <- glm(cbind(n_tb, n - n_tb) ~ com + hiv, 
                 family = binomial(), data = dat)
V_u_hiv   <- cluster_vcov(fit_u_hiv, dat)
terms_u_hiv <- grep("^hiv", names(coef(fit_u_hiv)), value = TRUE)
or_u_hiv  <- extract_or(fit_u_hiv, V_u_hiv, terms_u_hiv, hiv_levels[-1])
p_u_hiv   <- wald_test_terms(fit_u_hiv, V_u_hiv, terms_u_hiv)

# 13. Household TB Contact
fit_u_contact <- glm(cbind(n_tb, n - n_tb) ~ com + tb_contact, 
                     family = binomial(), data = dat)
V_u_contact   <- cluster_vcov(fit_u_contact, dat)
terms_u_contact <- grep("^tb_contact", names(coef(fit_u_contact)), value = TRUE)
or_u_contact  <- extract_or(fit_u_contact, V_u_contact, terms_u_contact, tb_contact_levels[-1])
p_u_contact   <- wald_test_terms(fit_u_contact, V_u_contact, terms_u_contact)

# 14. Residence and Mobility
fit_u_movement <- glm(cbind(n_tb, n - n_tb) ~ com + movement, 
                      family = binomial(), data = dat)
V_u_movement   <- cluster_vcov(fit_u_movement, dat)
terms_u_movement <- grep("^movement", names(coef(fit_u_movement)), value = TRUE)
or_u_movement  <- extract_or(fit_u_movement, V_u_movement, terms_u_movement, movement_levels[-1])
p_u_movement   <- wald_test_terms(fit_u_movement, V_u_movement, terms_u_movement)

#===============================================================================#
# SECTION 4: Fully Adjusted Model
#===============================================================================#
fit_adjusted <- glm(cbind(n_tb, n - n_tb) ~ com + agecat_2 + sex + hiv + tb_contact + movement,
                    family = binomial(),
                    data = dat)
fit_adjusted
V_adjusted <- cluster_vcov(fit_adjusted, dat)

#===============================================================================#
# SECTION 5: Extract Adjusted Results
#===============================================================================#
terms_a_com <- grep("^com", names(coef(fit_adjusted)), value = TRUE)
or_a_com  <- extract_or(fit_adjusted, V_adjusted, terms_a_com, community_levels[-1])
p_a_com   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_com)

terms_a_age <- grep("^agecat_2", names(coef(fit_adjusted)), value = TRUE)
or_a_age  <- extract_or(fit_adjusted, V_adjusted, terms_a_age, age_levels[-1])
p_a_age   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_age)

terms_a_sex <- grep("^sex", names(coef(fit_adjusted)), value = TRUE)
or_a_sex  <- extract_or(fit_adjusted, V_adjusted, terms_a_sex, sex_levels[-1])
p_a_sex   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_sex)

terms_a_hiv <- grep("^hiv", names(coef(fit_adjusted)), value = TRUE)
or_a_hiv  <- extract_or(fit_adjusted, V_adjusted, terms_a_hiv, hiv_levels[-1])
p_a_hiv   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_hiv)

terms_a_contact <- grep("^tb_contact", names(coef(fit_adjusted)), value = TRUE)
or_a_contact  <- extract_or(fit_adjusted, V_adjusted, terms_a_contact, tb_contact_levels[-1])
p_a_contact   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_contact)

terms_a_movement <- grep("^movement", names(coef(fit_adjusted)), value = TRUE)
or_a_movement  <- extract_or(fit_adjusted, V_adjusted, terms_a_movement, movement_levels[-1])
p_a_movement   <- wald_test_terms(fit_adjusted, V_adjusted, terms_a_movement)

#===============================================================================#
# SECTION 6: Construct Table Blocks
#===============================================================================#
make_table_block <- function(data, variable, characteristic, levels_order, unadjusted, adjusted, p_unadjusted, p_adjusted) {
  
  desc <- get_descriptive(data, variable, characteristic)
  
  reference <- tibble(level = levels_order[1],
                      OR_unadjusted = "Reference",
                      OR_adjusted   = "Reference")
  
  nonref <- tibble(level = levels_order[-1],
                   OR_unadjusted = sprintf("%.2f (%.2f, %.2f)", 
                                           unadjusted$estimate, 
                                           unadjusted$lower, 
                                           unadjusted$upper),
                   OR_adjusted   = sprintf("%.2f (%.2f, %.2f)", 
                                           adjusted$estimate, 
                                           adjusted$lower, 
                                           adjusted$upper))
  
  results <- bind_rows(reference, nonref)
  
  output <- desc %>%
    left_join(results, by = "level") %>%
    mutate(p_unadjusted = "", p_adjusted = "")
  
  output$p_unadjusted[1] <- if_else(p_unadjusted < 0.001, "<0.001", sprintf("%.3g", p_unadjusted))
  output$p_adjusted[1]   <- if_else(p_adjusted < 0.001, "<0.001", sprintf("%.3g", p_adjusted))
  
  output
}

block_com <- make_table_block(dat, "com", "Community", community_levels, or_u_com, or_a_com, p_u_com, p_a_com)
block_age <- make_table_block(dat, "agecat_2", "Age (in years)", age_levels, or_u_age, or_a_age, p_u_age, p_a_age)
block_sex <- make_table_block(dat, "sex", "Sex", sex_levels, or_u_sex, or_a_sex, p_u_sex, p_a_sex)
block_hiv <- make_table_block(dat, "hiv", "HIV/ART status", hiv_levels, or_u_hiv, or_a_hiv, p_u_hiv, p_a_hiv)
block_contact <- make_table_block(dat, "tb_contact", "Household TB contact", tb_contact_levels, or_u_contact, or_a_contact, p_u_contact, p_a_contact)
block_movement <- make_table_block(dat, "movement", "Residence and mobility status", movement_levels, or_u_movement, or_a_movement, p_u_movement, p_a_movement)

#===============================================================================#
# SECTION 7: Combine and Order Table
#===============================================================================#
table2_long <- bind_rows(block_com,
                         block_age,
                         block_sex,
                         block_hiv,
                         block_contact,
                         block_movement)

table2_long <- table2_long %>%
  mutate(
    characteristic_order = case_when(
      characteristic == "Community" ~ 1,
      characteristic == "Age (in years)" ~ 2,
      characteristic == "Sex" ~ 3,
      characteristic == "HIV/ART status" ~ 4,
      characteristic == "Household TB contact" ~ 5,
      characteristic == "Residence and mobility status" ~ 6
    ),
    
    level_order = case_when(
      characteristic == "Community" & level == "Z8"  ~ 1,
      characteristic == "Community" & level == "Z1"  ~ 2,
      characteristic == "Community" & level == "Z2"  ~ 3,
      characteristic == "Community" & level == "Z5"  ~ 4,
      characteristic == "Community" & level == "Z6"  ~ 5,
      characteristic == "Community" & level == "Z9"  ~ 6,
      characteristic == "Community" & level == "Z10" ~ 7,
      characteristic == "Community" & level == "Z11" ~ 8,
      
      characteristic == "Age (in years)" & level == "15 - 29" ~ 1,
      characteristic == "Age (in years)" & level == "30 - 39" ~ 2,
      characteristic == "Age (in years)" & level == "40 - 49" ~ 3,
      characteristic == "Age (in years)" & level == "50+"     ~ 4,
      
      characteristic == "Sex" & level == "Male"   ~ 1,
      characteristic == "Sex" & level == "Female" ~ 2,
      
      characteristic == "HIV/ART status" & level == "Tested negative by CHiP"                  ~ 1,
      characteristic == "HIV/ART status" & level == "Self-reported living with HIV & NOT on ART" ~ 2,
      characteristic == "HIV/ART status" & level == "Self-reported living with HIV & on ART"     ~ 3,
      characteristic == "HIV/ART status" & level == "Tested positive by CHiP"                  ~ 4,
      characteristic == "HIV/ART status" & level == "Not tested/indeterminate"                 ~ 5,
      
      characteristic == "Household TB contact" & level == "No"  ~ 1,
      characteristic == "Household TB contact" & level == "Yes" ~ 2,
      
      characteristic == "Residence and mobility status" & level == "Previously resident (same zone)" ~ 1,
      characteristic == "Residence and mobility status" & level == "Moved within community"         ~ 2,
      characteristic == "Residence and mobility status" & level == "Moved from outside community"    ~ 3
    )
  ) %>%
  arrange(characteristic_order, level_order)

table2 <- table2_long %>%
  mutate(
    `n/N (Yield)` = sprintf("%s/%s (%.1f)",
                            format(n_tb, big.mark = ",", trim = TRUE, scientific = FALSE),
                            format(N, big.mark = ",", trim = TRUE, scientific = FALSE),
                            yield)
  ) %>%
  select(characteristic, level, `n/N (Yield)`, OR_unadjusted, p_unadjusted, OR_adjusted, p_adjusted) %>%
  mutate(level = if_else(level == "50+", "≥50", level))

table2

#===============================================================================#
# SECTION 8: Create  Word Table
#===============================================================================#
table2_ft <- table2 %>%
  flextable() %>%
  set_header_labels(
    characteristic = "Potential Risk factors",
    level          = "",
    `n/N (Yield)`  = "n/N (Yield)",
    OR_unadjusted  = "OR (95%CI)",
    p_unadjusted   = "P-value*",
    OR_adjusted    = "AOR (95%CI)",
    p_adjusted     = "P-value*"
  ) %>%
  merge_v(j = "characteristic") %>%
  valign(j = "characteristic", valign = "top") %>%
  bold(j = "characteristic", bold = TRUE) %>%
  align(j = c("n/N (Yield)", "OR_unadjusted", "p_unadjusted", "OR_adjusted", "p_adjusted"), align = "center", part = "all") %>%
  add_header_row(
    values = c("", "", "Descriptive analysis", "Unadjusted analysis", "", "Adjusted analysis¶", ""),
    colwidths = c(1, 1, 1, 1, 1, 1, 1)
  ) %>%
  add_header_row(
    values = c("Round 3"),
    colwidths = 7
  ) %>%
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  autofit()

# Print Preview
print(table2_ft)

# Export Landscape Word Document
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape")
)

save_as_docx(
  table2_ft,
  path = "outputs/tables/Table_2_Associations.docx",
  pr_section = sect_properties
)


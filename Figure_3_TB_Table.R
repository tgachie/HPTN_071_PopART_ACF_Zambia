#===============================================================================#
# Script: 06_repeat_screening_analysis.R
# Purpose: Main R3 repeat-screening analysis and stable-resident sensitivity
# Inputs:  data/acf_repeat_screening.csv, data/acf_repeat_screening_stable.csv
# Outputs: outputs/tables/repeat_screening_main.{csv,docx}
#          outputs/tables/repeat_screening_sensitivity.{csv,docx}
#===============================================================================#
library(dplyr)
library(readr)
library(sandwich)
library(flextable)
library(officer)

if (!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)

#--- Helper Functions ----------------------------------------------------------#
cluster_vcov <- function(model, data) {
  V0 <- vcovCL(model, cluster = data$cluster_id, type = "HC0", cadjust = FALSE)
  N <- sum(data$n); G <- n_distinct(data$cluster_id); K <- length(coef(model))
  V0 * ((G / (G - 1)) * ((N - 1) / (N - K)))
}

wald_test_terms <- function(model, V, terms) {
  b <- coef(model)[terms]
  V_sub <- V[terms, terms, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V_sub) %*% b)
  pchisq(W, df = length(b), lower.tail = FALSE)
}

extract_or <- function(model, V, terms, levels_nonref) {
  beta <- coef(model)[terms]
  se <- sqrt(diag(V)[terms])
  tibble(level = levels_nonref, estimate = exp(beta),
         lower = exp(beta - 1.96 * se), upper = exp(beta + 1.96 * se))
}

format_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)

#--- Factor Levels -------------------------------------------------------------#
community_numeric_levels <- c(8, 1, 2, 5, 6, 9, 10, 11)
community_labels <- c("Z8", "Z1", "Z2", "Z5", "Z6", "Z9", "Z10", "Z11")
age_levels <- c("15 - 29", "30 - 39", "40 - 49", "50+")
sex_levels <- c("Male", "Female")
hiv_levels <- c("Self-reported living with HIV & on ART",
                "Self-reported living with HIV & NOT on ART",
                "Tested positive by CHiP", "Tested negative by CHiP",
                "Not tested/indeterminate")
tb_contact_levels <- c("No", "Yes")

#===============================================================================#
# MAIN ROUND 3 ANALYSIS
#===============================================================================#
main_dat <- read_csv("data/acf_repeat_screening.csv", show_col_types = FALSE)
repeat_levels_main <- c("First time screening", 
                        "Screened once before in the same zone", 
                        "Screened twice before in the same zone")

main_dat <- main_dat %>%
  mutate(
    repeat_screen = factor(repeat_screen, levels = repeat_levels_main),
    com = factor(com,levels = community_numeric_levels, labels = community_labels),
    agecat_2      = factor(agecat_2, levels = age_levels),
    sex           = factor(sex, levels = sex_levels),
    hiv           = factor(hiv, levels = hiv_levels),
    tb_contact    = factor(tb_contact, levels = tb_contact_levels),
    cluster_id    = interaction(com, zone, drop = TRUE)
  )
head(main_dat)

# Models
fit_main_1   <- glm(cbind(n_tb, n - n_tb) ~ repeat_screen + com, 
                    family = binomial(), data = main_dat)
V_main_1     <- cluster_vcov(fit_main_1, main_dat)
terms_main_1 <- grep("^repeat_screen", names(coef(fit_main_1)), value = TRUE)
or_main_1    <- extract_or(fit_main_1, V_main_1, terms_main_1, repeat_levels_main[-1])
p_main_1     <- wald_test_terms(fit_main_1, V_main_1, terms_main_1)

fit_main_2   <- glm(cbind(n_tb, n - n_tb) ~ repeat_screen + com + agecat_2 + sex + hiv + tb_contact, 
                    family = binomial(), data = main_dat)
V_main_2     <- cluster_vcov(fit_main_2, main_dat)
terms_main_2 <- grep("^repeat_screen", names(coef(fit_main_2)), value = TRUE)
or_main_2    <- extract_or(fit_main_2, V_main_2, terms_main_2, repeat_levels_main[-1])
p_main_2     <- wald_test_terms(fit_main_2, V_main_2, terms_main_2)

# Summaries & Formatting
main_desc <- main_dat %>%
  group_by(repeat_screen) %>%
  summarise(N = sum(n), n_tb = sum(n_tb), .groups = "drop") %>%
  mutate(yield = 100000 * n_tb / N, repeat_screen = as.character(repeat_screen))

main_reg1 <- bind_rows(tibble(level = repeat_levels_main[1], OR = "Reference"),
                       or_main_1 %>% transmute(level, OR = sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper)))
main_reg2 <- bind_rows(tibble(level = repeat_levels_main[1], AOR = "Reference"),
                       or_main_2 %>% transmute(level, AOR = sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper)))

main_results <- main_desc %>%
  rename(level = repeat_screen) %>%
  left_join(main_reg1, by = "level") %>%
  left_join(main_reg2, by = "level") %>%
  mutate(
    `n/N (Yield)`      = sprintf("%s/%s (%.1f)", format(n_tb, big.mark = ",", trim = TRUE), format(N, big.mark = ",", trim = TRUE), yield),
    `P-value`          = c(format_p(p_main_1), "", ""),
    `Adjusted P-value` = c(format_p(p_main_2), "", "")
  ) %>%
  select(`Prior screening history` = level, `n/N (Yield)`, OR, `P-value`, AOR, `Adjusted P-value`)

main_results

# Export Main
main_ft <- main_results %>%
  flextable() %>%
  set_header_labels(`Prior screening history` = "Prior screening history", `n/N (Yield)` = "n/N (Yield)", OR = "OR (95% CI)", `P-value` = "P-value", AOR = "AOR (95% CI)", `Adjusted P-value` = "P-value") %>%
  bold(part = "header") %>% align(align = "center", part = "header") %>% autofit()

sect_properties <- prop_section(page_size = page_size(orient = "landscape"))
save_as_docx(main_ft, path = "outputs/tables/repeat_screening_main.docx", pr_section = sect_properties)

#===============================================================================#
# STABLE-RESIDENT SENSITIVITY ANALYSIS
#===============================================================================#
stable_dat <- read_csv("data/acf_repeat_screening_stable.csv", show_col_types = FALSE)
repeat_levels_stable <- c("Screened once before in the same zone", 
                          "First time screening", 
                          "Screened twice before in the same zone")

community_numeric_levels <- c(8, 1, 2, 5, 6, 9, 10, 11)
community_labels <- c("Z8", "Z1", "Z2", "Z5", "Z6", "Z9", "Z10", "Z11")

#community_levels  <- c("Z8", "Z1", "Z2", "Z5", "Z6", "Z9", "Z10", "Z11")

stable_dat <- stable_dat %>%
  mutate(
    repeat_screen = factor(repeat_screen, levels = repeat_levels_stable),
    com = factor(com,levels = community_numeric_levels, labels = community_labels),
    #com           = factor(com, levels = community_levels),
    agecat_2      = factor(agecat_2, levels = age_levels),
    sex           = factor(sex, levels = sex_levels),
    hiv           = factor(hiv, levels = hiv_levels),
    tb_contact    = factor(tb_contact, levels = tb_contact_levels),
    cluster_id    = interaction(com, zone, drop = TRUE)
  )
head(stable_dat)
# Models
fit_stable_1   <- glm(cbind(n_tb, n - n_tb) ~ repeat_screen + com, 
                      family = binomial(), data = stable_dat)
V_stable_1     <- cluster_vcov(fit_stable_1, stable_dat)
terms_stable_1 <- grep("^repeat_screen", names(coef(fit_stable_1)), value = TRUE)
or_stable_1    <- extract_or(fit_stable_1, V_stable_1, terms_stable_1, repeat_levels_stable[-1])
p_stable_1     <- wald_test_terms(fit_stable_1, V_stable_1, terms_stable_1)

fit_stable_2   <- glm(cbind(n_tb, n - n_tb) ~ repeat_screen + com + agecat_2 + sex + hiv + tb_contact, 
                      family = binomial(), data = stable_dat)
V_stable_2     <- cluster_vcov(fit_stable_2, stable_dat)
terms_stable_2 <- grep("^repeat_screen", names(coef(fit_stable_2)), value = TRUE)
or_stable_2    <- extract_or(fit_stable_2, V_stable_2, terms_stable_2, repeat_levels_stable[-1])
p_stable_2     <- wald_test_terms(fit_stable_2, V_stable_2, terms_stable_2)

# Summaries & Formatting
stable_desc <- stable_dat %>%
  group_by(repeat_screen) %>%
  summarise(N = sum(n), n_tb = sum(n_tb), .groups = "drop") %>%
  mutate(yield = 100000 * n_tb / N, repeat_screen = as.character(repeat_screen))

stable_reg1 <- bind_rows(tibble(level = repeat_levels_stable[1], OR = "Reference"),
                         or_stable_1 %>% transmute(level, OR = sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper)))
stable_reg2 <- bind_rows(tibble(level = repeat_levels_stable[1], AOR = "Reference"),
                         or_stable_2 %>% transmute(level, AOR = sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper)))

stable_results <- stable_desc %>%
  rename(level = repeat_screen) %>%
  left_join(stable_reg1, by = "level") %>%
  left_join(stable_reg2, by = "level") %>%
  mutate(
    `n/N (Yield)`      = sprintf("%s/%s (%.1f)", format(n_tb, big.mark = ",", trim = TRUE), format(N, big.mark = ",", trim = TRUE), yield),
    `P-value`          = c(format_p(p_stable_1), "", ""),
    `Adjusted P-value` = c(format_p(p_stable_2), "", "")
  ) %>%
  select(`Prior screening history` = level, `n/N (Yield)`, OR, `P-value`, AOR, `Adjusted P-value`)
stable_results

# Export Sensitivity
stable_ft <- stable_results %>%
  flextable() %>%
  set_header_labels(`Prior screening history` = "Prior screening history", `n/N (Yield)` = "n/N (Yield)", OR = "OR (95% CI)", `P-value` = "P-value", AOR = "AOR (95% CI)", `Adjusted P-value` = "P-value") %>%
  bold(part = "header") %>% align(align = "center", part = "header") %>% autofit()

save_as_docx(stable_ft, path = "outputs/tables/repeat_screening_sensitivity.docx", pr_section = sect_properties)
    
    
    
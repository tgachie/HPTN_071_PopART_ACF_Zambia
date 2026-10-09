#===============================================================================#
# Script: 01_create_aggregated_data.R
# Purpose:
#   1. Import and clean the restricted individual-level PopART ACF dataset
#   2. Standardise variable coding and factor labels
#   3. Export non-identifiable aggregated datasets for open-access repository
#
# NOTE: The raw individual-level dataset ('data_private/') is restricted and 
#       MUST NOT be committed or pushed to public repositories.
#===============================================================================#

rm(list = ls())
#-------------------------------------------------------------------------------#
#1.  Libraries & Setup
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(readr)
library(labelled)

# Ensure local target directories exist
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

#-------------------------------------------------------------------------------#
# 2. Load restricted individual-level data
#-------------------------------------------------------------------------------#
dat <- read.csv("acf_minimal_dataset.csv")

#-------------------------------------------------------------------------------#
# 3. Standardise factor levels & set variable labels
#-------------------------------------------------------------------------------#
dat <- dat %>%
  mutate(
    ir = factor(ir, levels = 1:3, labels = c("Round 1", "Round 2", "Round 3")),
    sex = factor(sex, levels = 1:2, labels = c("Male", "Female")),
    agecat_1 = factor(agecat_1, levels = 1:8, 
                      labels = c("0 - 14", "15 - 17", "18 - 19", "20 - 24", 
                                 "25 - 29", "30 - 39", "40 - 49", "50+")),
    agecat_2 = factor(agecat_2, levels = 1:4, 
                      labels = c("15 - 29", "30 - 39", "40 - 49", "50+")),
    hiv = factor(hiv, levels = 1:5, 
                 labels = c("Self-reported living with HIV & NOT on ART",
                            "Tested positive by CHiP",
                            "Tested negative by CHiP",
                            "Not tested/indeterminate",
                            "Self-reported living with HIV & on ART")),
    movement = factor(movement, levels = 1:3, 
                      labels = c("Previously resident (same zone)",
                                 "Moved within community",
                                 "Moved from outside community")),
    seen = factor(seen, levels = 1:4, 
                  labels = c("In IR1 & IR3", "In IR1 to IR3", "In IR2 & IR3", "In IR3 only")),
    repeat_screen = factor(repeat_screen, levels = 1:3, 
                           labels = c("First time screening",
                                      "Screened once before in the same zone",
                                      "Screened twice before in the same zone")),
    across(c(tb_contact, enumerated, part_elig, participated, tb_treatment, 
             tb_screen, tb_screen_pos, sputum_collected, results_rec, tb, 
             started_TB_treatment),
           ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes")))
  ) %>%
  set_variable_labels(
    idno                 = "Participant unique ID across rounds 1 to 3",
    ir                   = "Annual round of the PopART intervention",
    com                  = "Community of residence",
    zone                 = "CHiP zone",
    sex                  = "Sex of participant",
    agecat_1             = "Age group, detailed",
    agecat_2             = "Age group used for TB yield analysis",
    hiv                  = "HIV/ART status",
    movement             = "Residence and mobility status",
    seen                 = "Residency/enumeration across intervention rounds",
    tb_contact           = "Household TB contact",
    enumerated           = "Enumerated",
    part_elig            = "Eligible to participate in intervention",
    participated         = "Participated in intervention",
    tb_treatment         = "Currently on TB treatment",
    tb_screen            = "Screened for TB symptoms",
    repeat_screen        = "Prior TB symptom-screening history during Round 3",
    tb_screen_pos        = "Screened positive for TB symptoms",
    sputum_collected     = "Sputum sample collected",
    results_rec          = "Diagnostic result recorded",
    tb                   = "Newly diagnosed TB",
    started_TB_treatment = "Started TB treatment"
  )

table(dat$tb_contact,dat$ir)
table(dat$movement,dat$ir)

#===============================================================================#
# SECTION 2: Create aggregated data for Figure 1 (Cascade)
#===============================================================================#
dat <- dat %>%
  mutate(
    below_age = case_when(
      ir == "Round 1" & agecat_1 %in% c("0 - 14", "15 - 17") ~ 1,
      ir %in% c("Round 2", "Round 3") & agecat_1 == "0 - 14" ~ 1,
      TRUE ~ 0
    )
  )

cascade_dat <- dat %>%
  group_by(ir, com, zone) %>%
  summarise(
    listed             = n(),
    below_age          = sum(below_age == 1, na.rm = TRUE),
    eligible           = sum(part_elig == "Yes", na.rm = TRUE),
    participated       = sum(participated == "Yes", na.rm = TRUE),
    on_tb_treatment    = sum(tb_treatment == "Yes", na.rm = TRUE),
    screened           = sum(tb_screen == "Yes", na.rm = TRUE),
    screen_positive    = sum(tb_screen_pos == "Yes", na.rm = TRUE),
    sputum_collected_n = sum(sputum_collected == "Yes", na.rm = TRUE),
    results_recorded   = sum(results_rec == "Yes", na.rm = TRUE),
    tb_detected        = sum(tb == "Yes", na.rm = TRUE),
    started_treatment  = sum(started_TB_treatment == "Yes", na.rm = TRUE),
    .groups            = "drop"
  )

write_csv(cascade_dat, "data/acf_cascade_by_zone.csv")

#===============================================================================#
# SECTION 3: Create aggregated data for Figure 2 (Subgroups)
#===============================================================================#
screened_dat <- dat %>% filter(tb_screen == "Yes")

# Helper function for Figure 2 aggregation
aggregate_fig2 <- function(data, group_var) {
  data %>%
    group_by(ir, com, zone, across(all_of(group_var))) %>%
    summarise(
      n_screened    = n(),
      n_symptom_pos = sum(tb_screen_pos == "Yes", na.rm = TRUE),
      n_tb          = sum(tb == "Yes", na.rm = TRUE),
      .groups       = "drop"
    )
}

fig2_sex <- aggregate_fig2(screened_dat, "sex")
fig2_age <- aggregate_fig2(screened_dat, "agecat_1")
fig2_hiv <- aggregate_fig2(screened_dat, "hiv")
fig2_tb_contact <- aggregate_fig2(screened_dat, "tb_contact")  #Not used in the figures but will be used in the Table_1_Summary

write_csv(fig2_sex, "data/acf_figure2_sex.csv")
write_csv(fig2_age, "data/acf_figure2_age.csv")
write_csv(fig2_hiv, "data/acf_figure2_hiv.csv")
write_csv(fig2_tb_contact, "data/fig2_tb_contact.csv")

#===============================================================================#
# SECTION 4: Create aggregated data for Regressions (Table 2)
#===============================================================================#
regression_r3 <- dat %>%
  filter(ir == "Round 3", tb_screen == "Yes") %>%
  group_by(com, zone, agecat_2, sex, hiv, tb_contact, movement) %>%
  summarise(
    n    = n(),
    n_tb = sum(tb == "Yes", na.rm = TRUE),
    .groups = "drop"
  )

write_csv(regression_r3, "data/acf_regression_r3.csv")

#===============================================================================#
# SECTION 5: Create aggregated dataset for Figure 3 (Cascade by Repeat Screening Status)
#===============================================================================#
# 5a. Round 3 repeat-screening analysis dataset (cascade)
figure3_dat <- dat %>%
  filter(
    ir %in% c("Round 3", 3),
    tb_screen %in% c("Yes", 1),
    !is.na(repeat_screen)
  ) %>%
  group_by(com, zone, repeat_screen) %>%
  summarise(
    screened = n(),
    screen_positive = sum(tb_screen_pos %in% c("Yes", 1), na.rm = TRUE),
    tb_detected = sum(tb %in% c("Yes", 1), na.rm = TRUE),
    started_treatment = sum(started_TB_treatment %in% c("Yes", 1), na.rm = TRUE),
    .groups = "drop"
  )

# Export to public aggregated data directory
write_csv(figure3_dat,"data/acf_figure3_cascade.csv")


# 5b. Round 3 repeat-screening regression analysis dataset
repeat_screening_dat <- dat %>%
  filter(ir == "Round 3", tb_screen == "Yes") %>%
  group_by(com, zone, agecat_2, sex, hiv, tb_contact, repeat_screen) %>%
  summarise(
    n    = n(),
    n_tb = sum(tb == "Yes", na.rm = TRUE),
    .groups = "drop"
  )

write_csv(repeat_screening_dat, "data/acf_repeat_screening.csv")


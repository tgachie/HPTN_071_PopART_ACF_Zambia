#===============================================================================#
# Script: 02_cascade_analysis.R
# Purpose: Replicate Figure 1 (TB Treatment Cascade Flowchart Summary)
# Input:   data/acf_cascade_by_zone.csv (Aggregated dataset)
# Output:  outputs/tables/cascade_table.docx
#          outputs/tables/cascade_table.csv
#===============================================================================#

#-------------------------------------------------------------------------------#
# 1. Setup & Libraries
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(readr)
library(flextable)
library(officer)

# Ensure output directory exists
if (!dir.exists("outputs/tables")) {
  dir.create("outputs/tables", recursive = TRUE)
}

#-------------------------------------------------------------------------------#
# 2. Load Aggregated Data
#-------------------------------------------------------------------------------#
data_path <- file.path("data", "acf_cascade_by_zone.csv")
cascade_raw <- read_csv(data_path, show_col_types = FALSE)

#-------------------------------------------------------------------------------#
# 3. Figure 1: Calculate Overall Cascade Counts & Percentages by Round
#-------------------------------------------------------------------------------#
cascade <- cascade_raw %>%
  group_by(ir) %>%
  summarise(
    listed             = sum(listed, na.rm = TRUE),
    below_age          = sum(below_age, na.rm = TRUE),
    eligible           = sum(eligible, na.rm = TRUE),
    participated       = sum(participated, na.rm = TRUE),
    on_tb_treatment    = sum(on_tb_treatment, na.rm = TRUE),
    screened           = sum(screened, na.rm = TRUE),
    screen_positive    = sum(screen_positive, na.rm = TRUE),
    sputum_collected_n = sum(sputum_collected_n, na.rm = TRUE),
    results_recorded   = sum(results_recorded, na.rm = TRUE),
    tb_detected        = sum(tb_detected, na.rm = TRUE),
    started_treatment  = sum(started_treatment, na.rm = TRUE),
    .groups            = "drop"
  ) %>%
  mutate(
    pct_below_age         = 100 * below_age / listed,
    pct_eligible          = 100 * eligible / listed,
    pct_participated      = 100 * participated / eligible,
    pct_on_tb_treatment   = 100 * on_tb_treatment / participated,
    pct_screened          = 100 * screened / participated,
    pct_screen_positive   = 100 * screen_positive / screened,
    pct_sputum_collected  = 100 * sputum_collected_n / screen_positive,
    pct_results_recorded  = 100 * results_recorded / sputum_collected_n,
    tb_yield_per_100000   = 100000 * tb_detected / screened,
    pct_started_treatment = 100 * started_treatment / tb_detected
  )

#-------------------------------------------------------------------------------#
# 4. Create Formatted Cascade Table
#-------------------------------------------------------------------------------#
cascade_table <- cascade %>%
  transmute(
    ir,
    `Household members listed`          = format(listed, big.mark = ","),
    `Below age eligibility threshold`  = sprintf("%s (%.1f%%)", format(below_age, big.mark = ","), pct_below_age),
    `Eligible to participate`          = sprintf("%s (%.1f%%)", format(eligible, big.mark = ","), pct_eligible),
    `Participated`                      = sprintf("%s (%.1f%%)", format(participated, big.mark = ","), pct_participated),
    `Already on TB treatment`           = sprintf("%s (%.2f%%)", format(on_tb_treatment, big.mark = ","), pct_on_tb_treatment),
    `Screened for TB symptoms`          = sprintf("%s (%.1f%%)", format(screened, big.mark = ","), pct_screened),
    `Screened positive for TB symptoms` = sprintf("%s (%.1f%%)", format(screen_positive, big.mark = ","), pct_screen_positive),
    `Sputum collected`                  = sprintf("%s (%.1f%%)", format(sputum_collected_n, big.mark = ","), pct_sputum_collected),
    `Results recorded`                  = sprintf("%s (%.1f%%)", format(results_recorded, big.mark = ","), pct_results_recorded),
    `MTB detected`                      = sprintf("%s (%.1f per 100,000 screened)", format(tb_detected, big.mark = ","), tb_yield_per_100000),
    `Started TB treatment`              = sprintf("%s (%.1f%%)", format(started_treatment, big.mark = ","), pct_started_treatment)
  ) %>%
  pivot_longer(cols = -ir, names_to = "Stage", values_to = "Value") %>%
  pivot_wider(names_from = ir, values_from = Value)

# Print console preview
print(cascade_table)

#-------------------------------------------------------------------------------#
# 5. Export Results
#-------------------------------------------------------------------------------#
# Define landscape page properties
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape")
)

# Export table as landscape Word document
cascade_table %>%
  flextable() %>%
  autofit() %>%
  save_as_docx(
    path = "outputs/tables/cascade_table.docx",
    pr_section = sect_properties
  )


#-------------------------------------------------------------------------------#
# 6. Figure 3: Calculate Overall Cascade Counts & Percentages by prior screening
#-------------------------------------------------------------------------------#








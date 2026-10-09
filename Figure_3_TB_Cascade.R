#===============================================================================#
# Script: 06_figure3_analysis.R
# Purpose: Generate Figure 3 Cascade Analysis by Prior Screening History
#
# Input:   data/acf_figure3_cascade.csv
#
# Outputs: outputs/tables/figure3_cascade_table.csv
#          outputs/tables/figure3_cascade_table.docx
#===============================================================================#

#-------------------------------------------------------------------------------#
# 1. Setup & Libraries
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(readr)
library(flextable)
library(officer)

# Ensure output directories exist
if (!dir.exists("outputs/tables")) {
  dir.create("outputs/tables", recursive = TRUE)
}

#-------------------------------------------------------------------------------#
# 2. Load Figure 3 Aggregated Data
#-------------------------------------------------------------------------------#
figure3_raw <- read_csv("data/acf_figure3_cascade.csv", show_col_types = FALSE)

# Restore screening-history factor ordering
repeat_screen_levels <- c(
  "First time screening", 
  "Screened once before in the same zone",
  "Screened twice before in the same zone"
)

figure3_raw <- figure3_raw %>%
  mutate(
    repeat_screen = factor(repeat_screen, levels = repeat_screen_levels)
  )

#-------------------------------------------------------------------------------#
# 3. Calculate Figure 3 Cascade Counts and Metrics
#-------------------------------------------------------------------------------#
figure3_cascade <- figure3_raw %>%
  group_by(repeat_screen) %>%
  summarise(
    screened          = sum(screened, na.rm = TRUE),
    screen_positive   = sum(screen_positive, na.rm = TRUE),
    tb_detected       = sum(tb_detected, na.rm = TRUE),
    started_treatment = sum(started_treatment, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  mutate(
    # Distribution of screening history among all screened in R3
    pct_screened = 100 * screened / sum(screened),
    
    # TB symptom prevalence within screening-history group
    pct_screen_positive = 100 * screen_positive / screened,
    
    # TB yield per 100,000 screened
    tb_yield_per_100000 = 100000 * tb_detected / screened,
    
    # Treatment initiation among TB cases in each screening-history group
    pct_started_treatment = 100 * started_treatment / tb_detected
  )

#-------------------------------------------------------------------------------#
# 4. Create Formatted Summary Table
#-------------------------------------------------------------------------------#
figure3_table <- figure3_cascade %>%
  transmute(
    `Prior screening history`  = as.character(repeat_screen),
    
    `Screened for TB symptoms` = sprintf(
      "%s (%.1f%%)",
      format(screened, big.mark = ",", scientific = FALSE, trim = TRUE),
      pct_screened
    ),
    
    `TB symptoms`              = sprintf(
      "%s (%.1f%%)",
      format(screen_positive, big.mark = ",", scientific = FALSE, trim = TRUE),
      pct_screen_positive
    ),
    
    `TB yield`                 = sprintf(
      "%s (%.1f per 100,000)",
      format(tb_detected, big.mark = ",", scientific = FALSE, trim = TRUE),
      tb_yield_per_100000
    ),
    
    `Started TB treatment`    = sprintf(
      "%s (%.1f%%)",
      format(started_treatment, big.mark = ",", scientific = FALSE, trim = TRUE),
      pct_started_treatment
    )
  )


#-------------------------------------------------------------------------------#
# 5. Create Flextable and Export to Landscape Word (.docx)
#-------------------------------------------------------------------------------#
figure3_ft <- figure3_table %>%
  flextable() %>%
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  align(j = c(2, 3, 4, 5), align = "center", part = "body") %>%
  autofit()

# Define landscape page properties
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape")
)

save_as_docx(
  figure3_ft,
  path = "outputs/tables/figure3_cascade_table.docx",
  pr_section = sect_properties
)


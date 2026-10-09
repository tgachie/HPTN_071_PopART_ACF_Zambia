#===============================================================================#
# Script: 03_figure2_analysis.R
# Purpose: Generate Figure 2 (8-panel plot of TB symptoms and yield)
# Inputs:  data/acf_figure2_sex.csv
#          data/acf_figure2_age.csv
#          data/acf_figure2_hiv.csv
# Output:  outputs/figures/TB_yield_Symptoms_panel_figure.tif
#===============================================================================#

#-------------------------------------------------------------------------------#
# 1. Setup & Libraries
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# Ensure output directory exists
if (!dir.exists("outputs/figures")) {
  dir.create("outputs/figures", recursive = TRUE)
}

#-------------------------------------------------------------------------------#
# 2. Helper Functions
#-------------------------------------------------------------------------------#
# Subpanel lettering helper
put.fig.letter <- function(label, location = "topleft", offset = c(0, 0), ...) {
  coords <- switch(
    location,
    topleft      = c(0.07, 0.97),
    topcenter    = c(0.5525, 0.98),
    topright     = c(0.985, 0.98),
    bottomleft   = c(0.015, 0.02),
    bottomcenter = c(0.5525, 0.02),
    bottomright  = c(0.985, 0.02),
    c(0.015, 0.98)
  )
  
  this.x <- grconvertX(coords[1] + offset[1], from = "nfc", to = "user")
  this.y <- grconvertY(coords[2] + offset[2], from = "nfc", to = "user")
  
  text(labels = label[1], x = this.x, y = this.y, xpd = TRUE, ...)
}

# Compute prevalence, yield, and exact 95% Binomial CIs
process_tb_data <- function(df) {
  df %>%
    mutate(
      id = row_number(),
      sympt_total = 100 * n_symptom_pos / n_screened,
      yield_total = 1e5 * n_tb / n_screened
    ) %>%
    mutate(
      ci_suspect = pmap(list(n_symptom_pos, n_screened), ~ binom.test(..1, ..2)$conf.int),
      ci_tb      = pmap(list(n_tb, n_screened), ~ binom.test(..1, ..2)$conf.int)
    ) %>%
    mutate(
      LCI_total_suspect = map_dbl(ci_suspect, 1) * 100,
      UCI_total_suspect = map_dbl(ci_suspect, 2) * 100,
      LCI_total_tb      = map_dbl(ci_tb, 1) * 1e5,
      UCI_total_tb      = map_dbl(ci_tb, 2) * 1e5
    ) %>%
    select(-ci_suspect, -ci_tb)
}

#-------------------------------------------------------------------------------#
# 3. Load & Process Aggregated Data
#-------------------------------------------------------------------------------#

df_sex <- read_csv("data/acf_figure2_sex.csv", show_col_types = FALSE)
df_age <- read_csv("data/acf_figure2_age.csv", show_col_types = FALSE)
df_hiv <- read_csv("data/acf_figure2_hiv.csv", show_col_types = FALSE)

# Map round levels to integer index (1, 2, 3)
map_round <- function(df) {
  df %>%
    mutate(
      round_num = case_when(
        ir == "Round 1" ~ 1,
        ir == "Round 2" ~ 2,
        ir == "Round 3" ~ 3,
        TRUE ~ as.numeric(ir)
      )
    )
}

df_sex <- map_round(df_sex)
df_age <- map_round(df_age)
df_hiv <- map_round(df_hiv)

# 3a. Overall Summary
df_overall <- df_sex %>%
  group_by(round_num) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  mutate(variable = "overall") %>%
  process_tb_data()

# 3b. Sex Summary
rs_male <- df_sex %>%
  filter(sex == "Male") %>%
  group_by(round_num) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  mutate(variable = "Male") %>%
  process_tb_data()

rs_female <- df_sex %>%
  filter(sex == "Female") %>%
  group_by(round_num) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  mutate(variable = "Female") %>%
  process_tb_data()

# 3c. Age Summary (Symptoms)
rs_age_symptom <- df_age %>%
  filter(!is.na(agecat_1)) %>%
  filter(!(round_num == 1 & agecat_1 == "15 - 17")) %>% # 15-17 excluded in Round 1
  group_by(round_num, agecat_1) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  rename(variable = agecat_1) %>%
  process_tb_data()

# 3d. Age Summary (Collapsed Yield: 15-29, 30-39, 40-49, 50+)
rs_age_yield <- df_age %>%
  filter(!is.na(agecat_1)) %>%
  filter(!(round_num == 1 & agecat_1 == "15 - 17")) %>%
  mutate(
    variable = if_else(agecat_1 %in% c("15 - 17", "18 - 19", "20 - 24", "25 - 29"), 
                       "15-29yrs", as.character(agecat_1))
  ) %>%
  mutate(variable = if_else(variable == "50+", "50+ yrs", variable)) %>%
  group_by(round_num, variable) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  process_tb_data()

# 3e. HIV Summary
rs_hiv <- df_hiv %>%
  filter(!is.na(hiv)) %>%
  group_by(round_num, hiv) %>%
  summarise(
    n_screened    = sum(n_screened, na.rm = TRUE),
    n_symptom_pos = sum(n_symptom_pos, na.rm = TRUE),
    n_tb          = sum(n_tb, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  rename(variable = hiv) %>%
  process_tb_data()

#-------------------------------------------------------------------------------#
# 4. Plotting Execution
#-------------------------------------------------------------------------------#

x_positions <- c(0.25, 1.00, 1.75)
x_labels    <- c("R1\n(2014–15)", "R2\n(2015–16)", "R3\n(2016–17)")

output_fig <- "outputs/figures/TB_yield_Symptoms_panel_figure.tif"

# Save figure to disk at 300 DPI
tiff(output_fig, width = 4000, height = 2000, res = 300)

par(mfrow = c(2, 4), mai = c(0.9, 0.75, 0.3, 0.35))

# PANEL (a): Overall TB symptoms
plot(x_positions, df_overall$sympt_total, ylim = c(0, 4), xlim = c(0, 2), pch = 15,
     ylab = "Prevalence of TB symptoms (%)", xaxt = "n", xlab = "",
     main = "Overall TB symptoms (by Round)")
arrows(x_positions, df_overall$LCI_total_suspect, x_positions, df_overall$UCI_total_suspect,
       length = 0.05, angle = 90, code = 3)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(a)", "topleft", font = 2, cex = 2)

# PANEL (b): TB symptoms by sex
male_x   <- x_positions
female_x <- x_positions + 0.1

plot(male_x, rs_male$sympt_total, ylim = c(0, 4), xlim = c(0, 2), pch = 15,
     ylab = "Prevalence of TB symptoms (%)", xaxt = "n", xlab = "",
     main = "TB symptoms Stratified by sex")
arrows(male_x, rs_male$LCI_total_suspect, male_x, rs_male$UCI_total_suspect, length = 0.05, angle = 90, code = 3)
points(female_x, rs_female$sympt_total, pch = 15, col = "blue")
arrows(female_x, rs_female$LCI_total_suspect, female_x, rs_female$UCI_total_suspect, col = "blue", length = 0.05, angle = 90, code = 3)
legend("topleft", legend = c("Male", "Female"), fill = c("black", "blue"), border = c("black", "blue"), cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(b)", "topleft", font = 2, cex = 2)

# PANEL (c): TB symptoms by age
age_vars    <- c("15 - 17", "18 - 19", "20 - 24", "25 - 29", "30 - 39", "40 - 49", "50+")
age_labels  <- c("15-17yrs", "18-19yrs", "20-24yrs", "25-29yrs", "30-39yrs", "40-49yrs", "50+ yrs")
age_colors  <- c("black", "red", "blue", "darkgreen", "purple", "orange", "brown")
age_offsets <- seq(-0.3, 0.3, length.out = length(age_vars))

plot(NA, xlim = c(0, 2.05), ylim = c(0, 5), xaxt = "n", xlab = "",
     ylab = "Prevalence of TB symptoms (%)", main = "TB symptoms Stratified by Age")
for (i in seq_along(age_vars)) {
  sub_data <- rs_age_symptom %>% filter(variable == age_vars[i]) %>% arrange(round_num)
  xvals <- x_positions[sub_data$round_num] + age_offsets[i]
  points(xvals, sub_data$sympt_total, pch = 15, col = age_colors[i])
  arrows(xvals, sub_data$LCI_total_suspect, xvals, sub_data$UCI_total_suspect,
         length = 0.05, angle = 90, code = 3, col = age_colors[i])
}
legend("topleft", legend = age_labels, fill = age_colors, border = age_colors, cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(c)", "topleft", font = 2, cex = 2)

# PANEL (d): TB symptoms by HIV status
hiv_vars   <- c("Tested positive by CHiP", 
                "Self-reported living with HIV & NOT on ART", 
                "Self-reported living with HIV & on ART", 
                "Tested negative by CHiP", 
                "Not tested/indeterminate")
hiv_labels <- c("Tested +ve", "SR +ve & not on ART", "SR +ve & on ART", "Tested -ve", "Not-tested/Indeterminate")
hiv_colors <- c("black", "red", "blue", "darkgreen", "brown")
offsets    <- seq(-0.2, 0.2, length.out = length(hiv_vars))

plot(NA, xlim = c(0, 2), ylim = c(0, 15), xaxt = "n", xlab = "",
     ylab = "Prevalence of TB symptoms (%)", main = "TB symptoms Stratified by HIV Status")
for (i in seq_along(hiv_vars)) {
  sub_data <- rs_hiv %>% filter(variable == hiv_vars[i]) %>% arrange(round_num)
  xvals <- x_positions[sub_data$round_num] + offsets[i]
  points(xvals, sub_data$sympt_total, pch = 15, col = hiv_colors[i])
  arrows(xvals, sub_data$LCI_total_suspect, xvals, sub_data$UCI_total_suspect,
         length = 0.05, angle = 90, code = 3, col = hiv_colors[i])
}
legend("topleft", legend = hiv_labels, fill = hiv_colors, border = hiv_colors, cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(d)", "topleft", font = 2, cex = 2)

# PANEL (e): Overall TB yield
plot(x_positions, df_overall$yield_total, ylim = c(0, 150), xlim = c(0, 2), pch = 15,
     ylab = "Yield of TB Disease/100 000", xaxt = "n", xlab = "",
     main = "Yield of TB Disease (by Round)")
arrows(x_positions, df_overall$LCI_total_tb, x_positions, df_overall$UCI_total_tb,
       length = 0.05, angle = 90, code = 3)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(e)", "topleft", font = 2, cex = 2)

# PANEL (f): TB yield by sex
plot(male_x, rs_male$yield_total, ylim = c(0, 250), xlim = c(0, 2), pch = 15,
     ylab = "Yield of TB Disease/100 000", xaxt = "n", xlab = "",
     main = "Yield of TB Disease Stratified by sex")
arrows(male_x, rs_male$LCI_total_tb, male_x, rs_male$UCI_total_tb, length = 0.05, angle = 90, code = 3)
points(female_x, rs_female$yield_total, pch = 15, col = "blue")
arrows(female_x, rs_female$LCI_total_tb, female_x, rs_female$UCI_total_tb, col = "blue", length = 0.05, angle = 90, code = 3)
legend("topleft", legend = c("Male", "Female"), fill = c("black", "blue"), border = c("black", "blue"), cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(f)", "topleft", font = 2, cex = 2)

# PANEL (g): TB yield by age
age_vars_yield    <- c("15-29yrs", "30 - 39", "40 - 49", "50+ yrs")
age_labels_yield  <- c("15-29yrs", "30-39yrs", "40-49yrs", "50+ yrs")
age_colors_yield  <- c("darkgreen", "purple", "orange", "brown")
age_yield_offsets <- seq(-0.2, 0.2, length.out = length(age_vars_yield))

plot(NA, xlim = c(0, 2.05), ylim = c(0, 350), xaxt = "n", xlab = "",
     ylab = "Yield of TB Disease/100 000", main = "Yield of TB Disease Stratified by Age")
for (i in seq_along(age_vars_yield)) {
  sub_data <- rs_age_yield %>% filter(variable == age_vars_yield[i]) %>% arrange(round_num)
  xvals <- x_positions[sub_data$round_num] + age_yield_offsets[i]
  points(xvals, sub_data$yield_total, pch = 15, col = age_colors_yield[i])
  arrows(xvals, sub_data$LCI_total_tb, xvals, sub_data$UCI_total_tb,
         length = 0.05, angle = 90, code = 3, col = age_colors_yield[i])
}
legend("topleft", legend = age_labels_yield, fill = age_colors_yield, border = age_colors_yield, cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(g)", "topleft", font = 2, cex = 2)

# PANEL (h): TB yield by HIV status
plot(NA, xlim = c(0, 2), ylim = c(0, 2000), xaxt = "n", xlab = "",
     ylab = "Yield of TB Disease/100 000", main = "Yield of TB Disease Stratified by HIV Status")
for (i in seq_along(hiv_vars)) {
  sub_data <- rs_hiv %>% filter(variable == hiv_vars[i]) %>% arrange(round_num)
  xvals <- x_positions[sub_data$round_num] + offsets[i]
  points(xvals, sub_data$yield_total, pch = 15, col = hiv_colors[i])
  arrows(xvals, sub_data$LCI_total_tb, xvals, sub_data$UCI_total_tb,
         length = 0.05, angle = 90, code = 3, col = hiv_colors[i])
}
legend("topleft", legend = hiv_labels, fill = hiv_colors, border = hiv_colors, cex = 0.8)
axis(1, at = x_positions, labels = FALSE)
mtext(x_labels, side = 1, at = x_positions, line = 2.5, cex = 0.9)
put.fig.letter("(h)", "topleft", font = 2, cex = 2)

dev.off()

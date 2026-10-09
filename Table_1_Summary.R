#===============================================================================#
# Script: 04_table1_analysis.R
# Purpose: Generate Table 1
#          Characteristics of participants screened for TB symptoms
#
# Inputs:
#   data/acf_figure2_sex.csv
#   data/acf_figure2_age.csv
#   data/acf_figure2_hiv.csv
#   data/acf_regression_r3.csv
#
# Outputs:
#   outputs/tables/Table_1_Summary.docx
#   outputs/tables/Table_1_Summary.csv
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
df_sex <- read_csv("data/acf_figure2_sex.csv", show_col_types = FALSE)
df_age <- read_csv("data/acf_figure2_age.csv", show_col_types = FALSE)
df_hiv <- read_csv("data/acf_figure2_hiv.csv", show_col_types = FALSE)
df_tbcontact <- read_csv("data/fig2_tb_contact.csv", show_col_types = FALSE)
df_r3  <- read_csv("data/acf_regression_r3.csv", show_col_types = FALSE)

#-------------------------------------------------------------------------------#
# 3. Restore Variable Levels and Ordering
#-------------------------------------------------------------------------------#
round_levels <- c("Round 1", "Round 2", "Round 3")
age_levels   <- c("15 - 17", "18 - 19", "20 - 24", "25 - 29", "30 - 39", "40 - 49", "50+")
sex_levels   <- c("Male", "Female")
hiv_levels   <- c(
  "Tested negative by CHiP",
  "Self-reported living with HIV & NOT on ART",
  "Self-reported living with HIV & on ART",
  "Tested positive by CHiP",
  "Not tested/indeterminate"
)
tb_contact_levels <- c("No", "Yes")
movement_levels   <- c(
  "Previously resident (same zone)",
  "Moved within community",
  "Moved from outside community"
)

# Apply factor ordering
df_age <- df_age %>%
  mutate(
    ir       = factor(ir, levels = round_levels),
    agecat_1 = factor(agecat_1, levels = age_levels)
  )

df_sex <- df_sex %>%
  mutate(
    ir  = factor(ir, levels = round_levels),
    sex = factor(sex, levels = sex_levels)
  )

df_hiv <- df_hiv %>%
  mutate(
    ir  = factor(ir, levels = round_levels),
    hiv = factor(hiv, levels = hiv_levels)
  )

df_tbcontact <- df_tbcontact %>%
  mutate(
    ir  = factor(ir, levels = round_levels),
    tb_contact = factor(tb_contact, levels = tb_contact_levels)
  )


df_r3 <- df_r3 %>%
  mutate(
    movement   = factor(movement, levels = movement_levels))

#-------------------------------------------------------------------------------#
# 4. Age Distribution Across Rounds
#-------------------------------------------------------------------------------#
t_age <- df_age %>%
  filter(!is.na(agecat_1)) %>%
  group_by(ir, agecat_1) %>%
  summarise(n = sum(n_screened, na.rm = TRUE), .groups = "drop") %>%
  group_by(ir) %>%
  mutate(
    denominator = sum(n),
    pct         = 100 * n / denominator
  ) %>%
  ungroup() %>%
  transmute(
    characteristic = "Age (in years)",
    level          = as.character(agecat_1),
    ir,
    n,
    denominator,
    pct
  )

#-------------------------------------------------------------------------------#
# 5. Sex Distribution Across Rounds
#-------------------------------------------------------------------------------#
t_sex <- df_sex %>%
  filter(!is.na(sex)) %>%
  group_by(ir, sex) %>%
  summarise(n = sum(n_screened, na.rm = TRUE), .groups = "drop") %>%
  group_by(ir) %>%
  mutate(
    denominator = sum(n),
    pct         = 100 * n / denominator
  ) %>%
  ungroup() %>%
  transmute(
    characteristic = "Sex",
    level          = as.character(sex),
    ir,
    n,
    denominator,
    pct
  )

#-------------------------------------------------------------------------------#
# 6. HIV/ART Status Across Rounds
#-------------------------------------------------------------------------------#
t_hiv <- df_hiv %>%
  filter(!is.na(hiv)) %>%
  group_by(ir, hiv) %>%
  summarise(n = sum(n_screened, na.rm = TRUE), .groups = "drop") %>%
  group_by(ir) %>%
  mutate(
    denominator = sum(n),
    pct         = 100 * n / denominator
  ) %>%
  ungroup() %>%
  transmute(
    characteristic = "HIV/ART status",
    level          = as.character(hiv),
    ir,
    n,
    denominator,
    pct
  )

#-------------------------------------------------------------------------------#
# 7. Household TB Contact Across Rounds
#-------------------------------------------------------------------------------#
t_contact <- df_tbcontact %>%
  filter(!is.na(tb_contact)) %>%
  group_by(ir, tb_contact) %>%
  summarise(n = sum(n_screened, na.rm = TRUE), .groups = "drop") %>%
  group_by(ir) %>%
  mutate(
    denominator = sum(n),
    pct         = 100 * n / denominator
  ) %>%
  ungroup() %>%
  transmute(
    characteristic = "Household TB contact",
    level          = as.character(tb_contact),
    ir,
    n,
    denominator,
    pct
  )

#-------------------------------------------------------------------------------#
# 8. Residence and Mobility Status - Round 3
#-------------------------------------------------------------------------------#
t_movement <- df_r3 %>%
  filter(!is.na(movement)) %>%
  group_by(movement) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    denominator = sum(n),
    pct         = 100 * n / denominator,
    ir          = factor("Round 3", levels = round_levels)
  ) %>%
  transmute(
    characteristic = "Residence and mobility status",
    level          = as.character(movement),
    ir,
    n,
    denominator,
    pct
  )

#-------------------------------------------------------------------------------#
# 9. Combine Table 1 Components
#-------------------------------------------------------------------------------#
table1_long <- bind_rows(t_age, t_sex, t_hiv, t_contact, t_movement)

#-------------------------------------------------------------------------------#
# 10. Define Variable and Category Ordering
#-------------------------------------------------------------------------------#
table1_long <- table1_long %>%
  mutate(
    # Variable order
    characteristic_order = case_when(
      characteristic == "Age (in years)"              ~ 1,
      characteristic == "Sex"                         ~ 2,
      characteristic == "HIV/ART status"              ~ 3,
      characteristic == "Household TB contact"        ~ 4,
      characteristic == "Residence and mobility status" ~ 5
    ),
    
    # Category order within each variable
    level_order = case_when(
      # Age
      characteristic == "Age (in years)" & level == "15 - 17" ~ 1,
      characteristic == "Age (in years)" & level == "18 - 19" ~ 2,
      characteristic == "Age (in years)" & level == "20 - 24" ~ 3,
      characteristic == "Age (in years)" & level == "25 - 29" ~ 4,
      characteristic == "Age (in years)" & level == "30 - 39" ~ 5,
      characteristic == "Age (in years)" & level == "40 - 49" ~ 6,
      characteristic == "Age (in years)" & level == "50+"     ~ 7,
      
      # Sex
      characteristic == "Sex" & level == "Male"   ~ 1,
      characteristic == "Sex" & level == "Female" ~ 2,
      
      # HIV/ART status
      characteristic == "HIV/ART status" & level == "Tested negative by CHiP"                  ~ 1,
      characteristic == "HIV/ART status" & level == "Self-reported living with HIV & NOT on ART" ~ 2,
      characteristic == "HIV/ART status" & level == "Self-reported living with HIV & on ART"     ~ 3,
      characteristic == "HIV/ART status" & level == "Tested positive by CHiP"                  ~ 4,
      characteristic == "HIV/ART status" & level == "Not tested/indeterminate"                 ~ 5,
      
      # Household TB contact
      characteristic == "Household TB contact" & level == "No"  ~ 1,
      characteristic == "Household TB contact" & level == "Yes" ~ 2,
      
      # Residence and mobility
      characteristic == "Residence and mobility status" & level == "Previously resident (same zone)" ~ 1,
      characteristic == "Residence and mobility status" & level == "Moved within community"         ~ 2,
      characteristic == "Residence and mobility status" & level == "Moved from outside community"    ~ 3
    ),
    
    # Retain explicit round ordering
    ir = factor(ir, levels = round_levels)
  ) %>%
  arrange(characteristic_order, level_order, ir)

#-------------------------------------------------------------------------------#
# 11. Format n (%)
#-------------------------------------------------------------------------------#
table1_long <- table1_long %>%
  mutate(
    value = sprintf(
      "%s (%.1f%%)",
      format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
      pct
    )
  )

#-------------------------------------------------------------------------------#
# 12. Convert to Wide Table and Preserve Round Order (Round 1, Round 2, Round 3)
#-------------------------------------------------------------------------------#
table1 <- table1_long %>%
  select(
    characteristic_order,
    level_order,
    characteristic,
    level,
    ir,
    value
  ) %>%
  pivot_wider(
    names_from = ir,
    values_from = value,
    names_expand = TRUE
  ) %>%
  arrange(characteristic_order, level_order) %>%
  select(
    characteristic,
    level,
    `Round 1`,
    `Round 2`,
    `Round 3`
  ) %>%
  mutate(
    across(
      c(`Round 1`, `Round 2`, `Round 3`),
      ~ replace_na(.x, "-")
    )
  )

#-------------------------------------------------------------------------------#
# 13. Calculate Total N for Round Column Headings
#-------------------------------------------------------------------------------#
round_n <- df_sex %>%
  group_by(ir) %>%
  summarise(N = sum(n_screened, na.rm = TRUE), .groups = "drop")

n_r1 <- round_n %>% filter(ir == "Round 1") %>% pull(N)
n_r2 <- round_n %>% filter(ir == "Round 2") %>% pull(N)
n_r3 <- round_n %>% filter(ir == "Round 3") %>% pull(N)

#-------------------------------------------------------------------------------#
# 15. Create Flextable
#-------------------------------------------------------------------------------#
table1_ft <- table1 %>%
  flextable() %>%
  set_header_labels(
    characteristic = "Participant characteristics",
    level          = "",
    `Round 1`      = paste0("Round 1\nN = ", format(n_r1, big.mark = ",", scientific = FALSE, trim = TRUE)),
    `Round 2`      = paste0("Round 2\nN = ", format(n_r2, big.mark = ",", scientific = FALSE, trim = TRUE)),
    `Round 3`      = paste0("Round 3\nN = ", format(n_r3, big.mark = ",", scientific = FALSE, trim = TRUE))
  ) %>%
  merge_v(j = "characteristic") %>%
  valign(j = "characteristic", valign = "top") %>%
  bold(j = "characteristic", bold = TRUE) %>%
  align(j = c("Round 1", "Round 2", "Round 3"), align = "left", part = "all") %>%
  autofit()

#-------------------------------------------------------------------------------#
# 16. Preview Table
#-------------------------------------------------------------------------------#
print(table1_ft)

#-------------------------------------------------------------------------------#
# 17. Export Table 1 as Landscape Word Document
#-------------------------------------------------------------------------------#
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape")
)

save_as_docx(
  table1_ft,
  path = "outputs/tables/Table_1_Summary.docx",
  pr_section = sect_properties
)

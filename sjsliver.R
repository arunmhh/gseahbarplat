library(readxl)
library(tidyr)
library(dplyr)

df <- read_excel("Rechenversion Laborparameter.xlsx", sheet = 1) %>%
  pivot_longer(cols = starts_with("dic"), names_to = "Attribute", values_to = "Value") %>%
  filter(Value != 0) %>%
  select(-c(starts_with("dic"))) %>%
  arrange(`BfArM Nummer`, Attribute)

unique_df <- unique(df[, c("BfArM Nummer", "Attribute", "Value")])
df_value <- df$Value
t.test(df_value)

# Download UBA db from https://www.umweltbundesamt.de/en/themen/chemikalien/arzneimittel/the-uba-database-pharmaceuticals-in-the-environment#data-aggregation
# Downloaded 20260510

library(readxl)
library(openxlsx)
library(tidyverse)

t.chemicals <- readRDS("t.chemicals.RDS")

raw.uba <- read_xlsx("data/pharms-uba_v3_2021_0.xlsx")
# correct headers are on row 2
col.head.uba <- raw.uba[2,]
clean.uba <- raw.uba %>% slice(3:n())
colnames(clean.uba) <- col.head.uba

# check how many we get get form this, only 91 cecs but tons of measreuments
sum(unique(clean.uba$`CAS number`) %in% t.chemicals$cas)

# extract the cols we need
uba.final <- clean.uba %>%
  filter(!is.na(`CAS number`), `CAS number` %in% t.chemicals$cas) %>%
  mutate(measurement_type = "MEC") %>%
  left_join(t.chemicals %>% select(cas, inchikey), join_by(`CAS number` == cas)) %>%
  select(
    inchikey,
    measurement_type,
    value_orig = `MEC original`,
    unit_orig = `Unit original`,
    value_stand = `MEC standardized`,
    unit_stand = `Unit standard`,
    country = `Sampling Country`,
    method = `Sampling Description`,
    timestamp_start = `Sampling Period Start`,
    timestamp_end = `Sampling Period End`,
    src = `Literature Citation`,
    src_desc = `Literature Type`
  ) %>%
  distinct()

write.xlsx(uba.final,
           "CEC_Table_Measurements_20260902.xlsx")

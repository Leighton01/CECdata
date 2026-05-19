# Download UBA db from https://www.umweltbundesamt.de/en/themen/chemikalien/arzneimittel/the-uba-database-pharmaceuticals-in-the-environment#data-aggregation
# Downloaded 20260510

library(readxl)
library(tidyverse)

raw.uba <- read_xlsx("data/pharms-uba_v3_2021_0.xlsx")
# correct headers are on row 2
col.head.uba <- raw.uba[2,]
clean.uba <- raw.uba %>% slice(3:n())
colnames(clean.uba) <- col.head.uba


sum(unique(clean.uba$`CAS number`) %in% merged.ordered$CAS)
# Only 16 are in the "suspected" list

# `CAS number`
# `Sampling Location`
# `MEC original`
# `Unit original`
# `MEC standardized`
# `Unit standard`
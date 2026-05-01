library(readxl)
library(tidyverse)
library(webchem)
library(janitor)



# Ping to see if api works
ping_service(
  service = c("bcpc", "chebi", "chembl", "cs", "cs_web", "cir", "cts", "etox", "fn",
              "nist", "opsin", "pc", "srs", "wd"),
  apikey = NULL
)

# Import and Clean --------------------------------------------------------
# Read working file (copy) retrieved on Marc 20, 2026
raw <- read_xlsx("data/Copy of CEC Database v1.2.xlsx")

# rename columns
raw.renamed <- raw
colnames(raw.renamed) <- raw[2,]
# remove last col (n/a) and first 2 rows (info)
raw.trimmed <- raw.renamed[1:(ncol(raw.renamed)-1)] %>% tail(., -2)
# Store the column names for eventual converion back
col.old <- colnames(raw.trimmed)
# rename crucial fields (CAS)
col.new <- col.old
col.new[col.new=="CAS Registry # (or EDF Substance ID)"] <- "CAS"
# replace col names
colnames(raw.trimmed) <- col.new
clean <- raw.trimmed
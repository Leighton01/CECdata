library(readxl)
library(tidyverse)
library(webchem)

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

# all that do not start with a number (ie EDF and... something else???)
clean %>% filter(!grepl("^[0-9]+", CAS))



# all that start with a number (ie CAS)
cas.val <- clean %>% filter(grepl("^[0-9]+", CAS))

# GET inchikey, but cas may correspond to multiple inchikey... :(
# OMG SO SLOW WHY
unlist(cts_convert(cas.val$CAS[1], from = "CAS", to = "InChIKey"))
inch <- unlist(cts_convert(cas.val$CAS, from = "CAS", to = "InChIKey",
            match = "first", verbose = FALSE))



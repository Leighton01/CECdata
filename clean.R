library(readxl)
library(tidyverse)
library(webchem)

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

# all that do not start with a number (ie EDF/others)
clean.edf <- clean %>% filter(!grepl("^[0-9]+", CAS))

# all that start with a number (ie CAS)
clean.cas <- clean %>% filter(grepl("^[0-9]+", CAS))

# Add inchkey to CAS
clean.cas$inch <- unlist(cts_convert(clean.cas$CAS, from = "CAS", to = "InChIKey",
            match = "first", verbose = FALSE))
# Need answer from email qs before doing the same to the EDF keys

# Get pubchem compound IDs based on inchikey, gotta join separately because it returns a tibble :(
cid <- get_cid(clean.cas$inch, from = "inchikey", match = "first", verbose = FALSE)
sum(is.na(cid$cid))
# clean.cas2 <- full_join(clean.cas, cid, by = c("inch" = "query"))


# Some inchikey not returning cid, though their CAS number does, might be some chemistry thing?
colnames(clean.cas2)
clean.cas2$cid
clean.cas2[1,]$inch
clean.cas2[35,]$inch

get_cid("33213-65-9", from = "inchikey", match = "first", verbose = FALSE)

# why not get pubchem cid directly from cas? why vignette suggests inchikey middle ground? check nas
clean.cas3 <- clean.cas
cid3 <- get_cid(clean.cas3$CAS, from = "cas", match = "first", verbose = FALSE)
# why are there NAs in the query/cas? there wasn't any in clean.cas3$CAS?
sum(is.na(clean.cas3$CAS))
sum(is.na(cid3$query))
filter(cid3, is.na(cid3$query))

temp <- left_join(clean.cas3, cid3, join_by("CAS"=="query"))
glimpse(temp)

# these are na because either it doesn't exist as CAS (some are pubchem SIDs) ,
# or the entry is actually a substance (thus not having a COMPOUND ID btu has a SID)
clean.cas4 <- (temp %>% filter(!is.na(cid)))
# check no na thank god
sum(is.na(clean.cas4$cid))



a <- pc_sect(1988, "Color / Form", domain = "compound", verbose=F)$Result


clean.cas5 <- clean.cas4


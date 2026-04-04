library(readxl)
library(tidyverse)
library(webchem)

# cid.map.cl2 is the reference df, 'Chemical Name' col + cid (retrieved),
#  no dupes or na

# clean.edf has all the cleaned entries that do not start with a number
# clean.cas has all the cleaned entries that do start with a number
# clean.cas2 is clean.cas + cid, filtered for rows with a valid cid
      # (ie no substances)

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

# Add identifiers (PubChem) ------------------------------------------------
## Inchikey Route (DEPRECATED) ---------------------------------------------------------
# # Add inchkey to CAS.. not sure if necessary? wait for email
# # Some inchikey not returning cid, though their CAS number does, might be some chemistry thing?
# clean.cas$inch <- unlist(cts_convert(clean.cas$CAS, from = "CAS", to = "InChIKey",
#             match = "first", verbose = FALSE))
# # Need answer from email qs before doing the same to the EDF keys


## CAS Route (DEPRECATED) --------------------------------------------------------------
# Split into CAS and others
# all that do not start with a number (mostly EDF)
clean.edf <- clean %>% filter(!grepl("^[0-9]+", CAS))

# all that start with a number (mostly CAS)
clean.cas <- clean %>% filter(grepl("^[0-9]+", CAS))

# get pubchem cid directly from cas
# cid <- get_cid(clean.cas$CAS, from = "cas", domain = "compound",
#                match = "first", verbose = FALSE)

# When there are multiple CIDs, get one with highest LiteratureCount
# not perfect, but good enough?
# cid.map <- get_cid(
#   clean.cas$CAS,
#   from = "cas",
#   domain = "compound",
#   match = "all",
#   verbose = FALSE
# )
#
# #
# props <- pc_prop(
#   cid_map$cid,
#   properties = "LiteratureCount"
# )
#
# saveRDS(cid.map, file="cid.map.rds")
# saveRDS(props, file="props.rds")
# load("cid.map.rds")
# load("props.rds")
#
#

# props.max <- props %>% group_by(CID) %>%
#   mutate(LiteratureCount = coalesce(LiteratureCount, -Inf)) %>%
#   slice_max(order_by = LiteratureCount, n = 1)
#
# cid.best <- cid.map %>%
#   left_join(props.max, join_by("cid"=="CID")) %>%
#   ungroup()

# # Some values return na for CID because it's not CAS (actually pubchem SIDs) ,
# # or the CAS is actually for a substance (thus not having a COMPOUND ID but SID)

## Name Route -------------------------------------------------------------
# Obtain CID from chemical name col (may contain name or formula, but both work)
cid.map2 <- get_cid(
  clean$`Chemical Name`,
  match = "all",
  verbose = FALSE)

saveRDS(cid.map2, file="cid.map2.rds")

# surprisingly, works better than CAS, returns only 1 CID for all compounds
# substances and ill-formatted names return NA (69 NAs)
cid.map2 %>% filter(is.na(cid)) %>% length

cid.map.cl <- cid.map2 %>% filter(!is.na(cid))
saveRDS(cid.map.cl, file="cid.map.cl.rds")
# There are some duplicates, such as ALPHA&BETA, CIS&TRANS, BETA&GAMMA,
# yields same result when searached manually? not yet filled
# Only 3 pairs, will remove for now

library(janitor)
dupes <- get_dupes(cid.map.cl, cid)$cid
cid.map.cl2 <- cid.map.cl %>% filter(!(cid %in% dupes))
saveRDS(cid.map.cl2, file="cid.map.cl2.rds")






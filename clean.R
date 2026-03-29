library(readxl)
library(tidyverse)
library(webchem)

# clean.edf has all the cleaned entries that does not start with a number
# clean.cas has all the cleaned entries that does start with a number
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

# Split into CAS and others ------------------------------------------
# all that do not start with a number (mostly EDF)
clean.edf <- clean %>% filter(!grepl("^[0-9]+", CAS))

# all that start with a number (mostly CAS)
clean.cas <- clean %>% filter(grepl("^[0-9]+", CAS))


# Add identifiers (PubChem) ------------------------------------------------
# # Add inchkey to CAS.. not sure if necessary? wait for email
# # Some inchikey not returning cid, though their CAS number does, might be some chemistry thing?
# clean.cas$inch <- unlist(cts_convert(clean.cas$CAS, from = "CAS", to = "InChIKey",
#             match = "first", verbose = FALSE))
# # Need answer from email qs before doing the same to the EDF keys

# get pubchem cid directly from cas
# cid <- get_cid(clean.cas$CAS, from = "cas", domain = "compound",
#                match = "first", verbose = FALSE)

# When there are multiple CIDs, get one with highest  LiteratureCount
# not perfect, but good enough?
cid_map <- get_cid(
  clean.cas$CAS,
  from = "cas",
  domain = "compound",
  match = "all",
  verbose = FALSE
)

#
props <- pc_prop(
  cid_map$cid,
  properties = "LiteratureCount"
)

saveRDS(cid_map, file="cid_map.rds")
saveRDS(props, file="props.rds")

cid_best <- cid_map %>%
  left_join(props, join_by("cid"=="CID")) %>%
  mutate(LiteratureCount = coalesce(LiteratureCount, -Inf)) %>%
  group_by(query) %>%
  slice_max(LiteratureCount, n = 1, with_ties = FALSE) %>%
  ungroup()



# Some values return na for CID because it's not CAS (actually pubchem SIDs) ,
# or the CAS is actually for a substance (thus not having a COMPOUND ID but SID)

# Add and filter for records with valid CIDs
clean.cas2 <- left_join(clean.cas, cid, join_by("CAS"=="query")) %>%
                  filter(!is.na(cid))


# Scrape for other properties (PubChem) -----------------------------------
## Description -------------------------------------------------------------
# Agrochemical Information or Biologic Information? try
a1 <- pc_sect(1988, "Agrochemical Information", domain = "compound", verbose=F)
a2 <- pc_sect(1988, "Biologic Information", domain = "compound", verbose=F)

## Physical Properties -----------------------------------------------------
# look at return and concatenate, connection often timeout
# # try with 1; success
# a1 <- pc_sect(1988, "Color / Form", domain = "compound", verbose=F)
# a2 <- pc_sect(1988, "Odor", domain = "compound", verbose=F)
# a3 <- pc_sect(1988, "Boiling Point", domain = "compound", verbose=F)
# paste(unlist(list(a1$Result, a2$Result, paste("Boiling Point: ", a3$Result))), collapse = ". ")
#
# # try with 2; success
# b1 <- pc_sect(c(1988,176), "Color / Form", domain = "compound", verbose=F)
# b2 <- pc_sect(c(1988,176), "Odor", domain = "compound", verbose=F)
# b3 <- pc_sect(c(1988,176), "Boiling Point", domain = "compound", verbose=F)
#
# b0 <- bind_rows(b1,b2,b3)
# b <- b0 %>% group_by(CID) %>% summarise(paste(Result, collapse=". "))

pp.colour <- pc_sect(clean.cas2$cid, "Color / Form", domain = "compound", verbose=F)

pp.odor <- pc_sect(clean.cas2$cid, "Odor", domain = "compound", verbose=F)
pp.bp <- pc_sect(clean.cas2$cid, "Boiling Point", domain = "compound", verbose=F)
pp.mp <- pc_sect(clean.cas2$cid, "Melting Point", domain = "compound", verbose=F)
pp.decompn <- pc_sect(clean.cas2$cid, "Decomposition", domain = "compound", verbose=F)

# Not sure why some are NA? sampled a few and they do return a result??
# pc_sect(2078, "Decomposition", domain = "compound", verbose=F)
# try as numeric ?
c(
  sum(is.na(pp.colour$Result)),
  sum(is.na(pp.odor$Result)),
  sum(is.na(pp.bp$Result)),
  sum(is.na(pp.mp$Result)),
  sum(is.na(pp.decompn$Result))
)

na.colour.cid
pp.colour

# the specific cids from cas does not return anything even though...
# some have other cids that MIGHT have the info, last match better?
get_cid("116-06-3", from = "cas", domain="compound", match = "all", verbose = FALSE)
pc_sect(3224, "Color / Form", verbose=F)

try <- get_cid("21259-20-1", from = "cas", domain="compound", match = "all", verbose = FALSE)

try.prop <- pc_prop(unlist(try$cid))

try.prop %>% select()
cid.best <- try.prop$CID[
  which.max(replace(try.prop$LiteratureCount, is.na(try.prop$LiteratureCount), -Inf))
]
# rank by LiteratureCount?

#
# try.syn <- pc_synonyms(unlist(try$cid), from = "cid", match = "all", verbose = F)
# # or go by most synonyms lol


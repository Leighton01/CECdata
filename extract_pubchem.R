library(readxl)
library(tidyverse)
library(webchem)

# Ping to see if api works
ping_service(
  service = c("bcpc", "chebi", "chembl", "cs", "cs_web", "cir", "cts", "etox", "fn",
              "nist", "opsin", "pc", "srs", "wd"),
  apikey = NULL
)

## Get CID from CAS --------------------------------------------------------------
# cid.map.cl2 is the reference df, 'Chemical Name' col + cid (retrieved),
#  no dupes or na

# clean.edf has all the cleaned entries that do not start with a number
# clean.cas has all the cleaned entries that do start with a number
# clean.cas2 is clean.cas + cid, filtered for rows with a valid cid
# (ie no substances)
# Split into CAS and others
# all that do not start with a number (mostly EDF)
clean.edf <- clean %>% filter(!grepl("^[0-9]+", CAS))

# all that start with a number (mostly CAS)
clean.cas <- clean %>% filter(grepl("^[0-9]+", CAS))

# get pubchem cid directly from cas
cas.to.cid <- get_cid(clean.cas$CAS, from = "cas", domain = "compound",
               match = "first", verbose = FALSE)



# # When there are multiple CIDs, get one with highest LiteratureCount
# # not perfect, but good enough?
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


# props.max <- props %>% group_by(CID) %>%
#   mutate(LiteratureCount = coalesce(LiteratureCount, -Inf)) %>%
#   slice_max(order_by = LiteratureCount, n = 1)
#
# cid.best <- cid.map %>%
#   left_join(props.max, join_by("cid"=="CID")) %>%
#   ungroup()

# # Some values return na for CID because it's not CAS (actually pubchem SIDs) ,
# # or the CAS is actually for a substance (thus not having a COMPOUND ID but SID)

## Get CID from Chemical Name-------------------------------------------------
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

# dupes need manual attention
dupes <- get_dupes(cid.map.cl, cid)
cid.map.cl2 <- cid.map.cl %>% filter(!(cid %in% dupes))

saveRDS(cid.map.cl2, file="cid.map.cl2.rds")

# Remove the chemicals that are in the dupe list (those need manual attnetion)
diff.pc.name1 <- setdiff(clean$`Chemical Name`, merged.ordered$`Chemical Name`)
diff.pc.name2 <- diff.pc.name1[which(!(diff.pc.name1 %in% dupes$query))]

# Add CAS to the list
diff.pc1 <- tibble(`Chemical Name`=diff.pc.name2)
diff.pc2 <- diff.pc1 %>% left_join(clean %>% select(`Chemical Name`, CAS),
                                   join_by(`Chemical Name` == `Chemical Name`))

# Remove anything that's don't start with a number and contains a dash
# (format of CAS)
diff.pc3 <- diff.pc2 %>% filter(grepl("^[0-9]+.*-", CAS))

# Look up inchikey from this list, using CAS
cas.inchi1 <- cts_convert(diff.pc3$CAS, from = "CAS", to = "PubChem CID")
cas.inchi2 <- cts_convert(diff.pc3$CAS, from = "CAS", to = "InChIKey")



# Retrieve Column Data ----------------------------------------------------
# ID
id_all_type <- pc_prop(cid.map.cl2$cid, properties = c("InChIKey","InChI",
                                                       "SMILES","IUPACName",
                                                       "MolecularFormula"
                                                       ), verbose = F)
# Chemical Properties

desc2 <- pc_sect(cid.map.cl2$cid, "Record Description", domain = "compound", verbose=F)

# Physical Properties
pp.colour <- pc_sect(cid.map.cl2$cid, "Color / Form", domain = "compound", verbose=F)
pp.odor <- pc_sect(cid.map.cl2$cid, "Odor", domain = "compound", verbose=F)
pp.bp <- pc_sect(cid.map.cl2$cid, "Boiling Point", domain = "compound", verbose=F)
pp.mp <- pc_sect(cid.map.cl2$cid, "Melting Point", domain = "compound", verbose=F)
pp.decompn <- pc_sect(cid.map.cl2$cid, "Decomposition", domain = "compound", verbose=F)


# some compounds don't have all sections, that's ok
c(
  sum(is.na(pp.colour$Result)),
  sum(is.na(pp.odor$Result)),
  sum(is.na(pp.bp$Result)),
  sum(is.na(pp.mp$Result)),
  sum(is.na(pp.decompn$Result))
)

# Remove duplicates
pp.colour2 <- pp.colour %>% distinct(CID, Result) %>% filter(!is.na(Result))
pp.odor2 <- pp.odor %>% distinct(CID, Result) %>% filter(!is.na(Result))
pp.bp2 <- pp.bp %>% distinct(CID, Result) %>% filter(!is.na(Result))
pp.mp2 <- pp.mp %>% distinct(CID, Result) %>% filter(!is.na(Result))
pp.decompn2 <- pp.decompn %>% distinct(CID, Result) %>% filter(!is.na(Result))

# Collapse and add tag
pp.bp3 <- pp.bp2 %>% group_by(CID) %>%
  summarise(Result= paste(Result, collapse = "; "), .groups = "drop")
pp.bp3$Result <- paste("Reported BPs: ", pp.bp3$Result)

pp.mp3 <- pp.mp2 %>% group_by(CID) %>%
  summarise(Result= paste(Result, collapse = "; "), .groups = "drop")
pp.mp3$Result <- paste("Reported MPs: ", pp.mp3$Result)

# Combine all pp dfs
pp <- rbind(pp.colour2, pp.odor2, pp.bp3, pp.mp3, pp.decompn2)

pp.comb <- pp %>% group_by(CID) %>%
  summarise(Result = paste(Result, collapse = ".\n"), .groups = "drop")

# saveRDS(pp.comb, "pp.comb.RDS")
#
# pp.by.name <- pp.comb %>% left_join(cid.map.cl, join_by (CID == cid)) %>% glimpse
# write.csv(pp.by.name, "pp.name.csv", fileEncoding="Windows-1252", row.names = FALSE)



# Chemical Releases

fate.exp <- pc_sect(cid.map.cl2$cid, "Environmental Fate/Exposure Summary",
                domain = "compound", verbose=F)

# Environmental Effects

fate <- pc_sect(cid.map.cl2$cid, "Environmental Fate",
                    domain = "compound", verbose=F)

# Human Health Effects
signs <- pc_sect(cid.map.cl2$cid, "Signs and Symptoms",
                 domain = "compound", verbose=F)
# group.save(signs, "signs.RDS")

# Sources
sources <- pc_sect(cid.map.cl2$cid, "Sources/Uses",
                   domain = "compound", verbose=F)

# Chemical Use (Source)
# Sources Facility / Location
# Source Industry

# Drinking Water Sources and Watersheds
water <- pc_sect(cid.map.cl2$cid, "Environmental Water Concentrations",
        domain = "compound", verbose=F)
# group.save(water)

# Drinking Water Sources
# drink <-

# Receiving Watersheds, paragraph 2
# shed <-


# Chemical Entry Points

# 12.2.8 Environmental Fate (S2P1)


# Monitoring Requirements
req <- pc_sect(cid.map.cl2$cid, "Regulatory Information",
                 domain = "compound", verbose=F)

# Removal Technologies
tech <- pc_sect(cid.map.cl2$cid, "Environmental Biodegradation",
               domain = "compound", verbose=F)

# Safe Production
safe <- pc_sect(cid.map.cl2$cid, "Storage Conditions",
                domain = "compound", verbose=F)
# Safe Use
uses <- pc_sect(cid.map.cl2$cid, "Personal Protective Equipment",
               domain = "compound", verbose=F)

# Safe Disposal
disposal <- pc_sect(cid.map.cl2$cid, "Disposal Methods",
               domain = "compound", verbose=F)

# Consumer Products
consum <- pc_sect(cid.map.cl2$cid, "Household Products",
                    domain = "compound", verbose=F)

# Exposure Routes

exp.routes <-  pc_sect(cid.map.cl2$cid, "Exposure Routes",
                          domain = "compound", verbose=F)

# Exposure Baseline

exp.base <- pc_sect(cid.map.cl2$cid, "Metabolism / Metabolites",
                          domain = "compound", verbose=F)
# 8.3 Metabolism / Metabolites (P2)


# Transgenerational Effects

generation <- pc_sect(cid.map.cl2$cid, "Health Effects",
                        domain = "compound", verbose=F)

# Hormetic Effects

carcin <- pc_sect(cid.map.cl2$cid, "Evidence for Carcinogenicity",
                      domain = "compound", verbose=F)

# H2O Sol

h20 <- pc_sect(cid.map.cl2$cid, "Volatilization from Water / Soil",
                      domain = "compound", verbose=F)



# Combine -----------------------------------------------------------------
list1 <- lst( "Chemical Properties" = desc2,
              "Physical Properties" = pp.comb,
              "Chemical Releases" = fate.exp,
              "Environmental Effects" = fate,
              "Human Health Effects" = signs,
              # water, not sure how to divide this data yet
              "Monitoring Requirements" = req,
              "Removal Technologies" = tech,
              "Safe Production" = safe,
              # "Safe Use" = uses, # no result
              #"Safe Disposal" = disposal, not sure how to divide this data yet
              "Consumer Products" = consum,
              "Exposure Routes" = exp.routes,
              #"Exposure Baseline" = exp.base, not sure how to divide this data yet
              "Transgenerational Effects" = generation,
              "Hormetic Effects" = carcin,
              "H2O Sol." = h20)

# check if any is empty
which((vapply(list1, function(x) !("Result" %in% colnames(x)), logical(1))))

# Drop unneeded columns
list2 <- lst()

# Clean
for (i in seq_along(list1)){

  # if (ncol(list1[[i]]) == 0) {
  #   temp_cn <- names(list1)[[i]]
  #   list2[[i]] <- list1[[i]]
  #   next
  # }

#  remove unneeded cols
 df.comb1 <- list1[[i]] %>% select("CID", "Result")

# drop na, group content by CID
 df.comb2 <- df.comb1 %>%
   tidyr::drop_na(Result) %>% group_by(CID) %>%
   summarise(Result = paste(Result, collapse = ".\n"), .groups = "drop")

 #  rename based on tbl names
 colnames(df.comb2) <- c("CID", names(list1)[i])
 list2[[i]] <- df.comb2

}

# check the resulting dfs
sapply(list2, function(x) nrow(x))
lapply(list2, function(x) colnames(x))

saveRDS(list2, "list2.RDS")
# load("list2.RDS")

# Merge
merged.pc.df <- purrr::reduce(list2, full_join, by = "CID")
# Check all cids are unique
sum(duplicated(merged.pc.df))

# Get CAS from orig file
cas.all <- clean %>% select(`Chemical Name`, CAS)

# Add CAS to cid map
cid.cas.map <- left_join(cid.map.cl, cas.all, join_by(query == 'Chemical Name'))

# Add chemical name + CAS to db
merged.pc.name <- merged.pc.df %>% left_join(cid.cas.map, join_by (CID == cid))
colnames(merged.pc.name)[which(colnames(merged.pc.name) == "query")] <- "Chemical Name"

# Add other IDs
merged.pc.id <- merged.pc.name %>% left_join(id_all_type, join_by (CID == CID))

# Reorder
merged.ordered <- merged.pc.id %>%
  select(`Chemical Name`, InChIKey, CAS, CID,
          CID, IUPACName, InChI, SMILES,
          MolecularFormula, everything()) %>%
  arrange(merged.ordered, `Chemical Name`)


# Save
openxlsx::write.xlsx(merged.ordered, "merged.ordered.xlsx")
# writexl::write_xlsx(merged.pc.name, "write.xlsx")



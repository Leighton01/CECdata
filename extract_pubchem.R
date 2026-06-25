library(readxl)
library(openxlsx)
library(tidyverse)
library(webchem)
library(janitor)
library(ctxR)


# CTS not working since beginning of May 2026 due to recent cyber attack -May 17 2026
# run cts relevant lines once fixed

# id.all2 contains records that have either CID or inchikey or cas (cid.map.cl2 has CID from name, cas.id has CID or inchikey from CAS)
# attn contains records that require manual attention

# Ping to see if api works
ping_service(
  service = c("bcpc", "chebi", "chembl", "cs", "cs_web", "cir", "cts", "etox", "fn",
              "nist", "opsin", "pc", "srs", "wd"),
  apikey = NULL
)

# check chemspider api key
cs_check_key()

# register epa api key, store in envir
register_ctx_api_key("4d4cbbce-4265-4b24-9712-326c0d8ea618", write = T)

# TABLE: CHEMICALS --------------------------------------------------------
# # Get InChIKey from CAS, CIR api
# inch.map.cas.cir <- cir_query(clean$CAS,
#                                representation = "stdinchikey",
#                                match = "all",
#                                verbose = FALSE)
#
# saveRDS(inch.map.cas.cir, "inch.map.cas.cir_0608.RDS")
inch.map.cas.cir <- readRDS("inch.map.cas.cir_0608.RDS")

# Many CAS will have multiple inchikeys. Add only uniques to list.
# var naming structure: "inchikey from cas, api CIR, unique only"
inch.map.cas.cir.uniq <- inch.map.cas.cir %>%
  # group by cas
  group_by(cas = query) %>%
  # keep if only 1 record and not na
  filter(n()==1, !(is.na(stdinchikey))) %>%
  ungroup() %>%
  rename(inchikey = stdinchikey) %>%
  # clean the cir inchikey results
  mutate(inchikey = str_replace(inchikey, "InChIKey=", ""))

# combine clean + inchikeys
clean.inch <- inch.map.cas.cir.uniq %>%
  left_join(clean,
            join_by(query == cas))%>%
  reframe(name, cas, inchikey, label)

# # Get CID from inchikey (retrieved through CIR) from pubchem API
# cid.map.inch.pc <- get_cid(
#   inch.map.cas.cir.uniq$InChIKey,
#   from = "inchikey",
#   match = "all",
#   verbose = FALSE)
#
# saveRDS(cid.map.inch.pc, "cid.map.inch.pc_0608.RDS")
cid.map.inch.pc <- readRDS("cid.map.inch.pc_0608.RDS")

# Again, keep only the unique CIDs
cid.map.inch.pc.uniq <- cid.map.inch.pc %>%
  group_by(query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  rename(inchikey = query) %>%
  left_join(clean.inch, join_by(inchikey == inchikey))

# # Get CID from chemical name
# cid.map.name.pc <- get_cid(
#   clean$name,
#   from = "name",
#   match = "all",
#   verbose = FALSE)
#
# saveRDS(cid.map.name.pc, "cid.map.name.pc_0608.RDS")
cid.map.name.pc <- readRDS("cid.map.name.pc_0608.RDS")

cid.map.name.pc.uniq <- cid.map.name.pc %>%
  group_by(query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  rename(name = query) %>%
  left_join(clean.inch, join_by(name == name))


# # Get CID from CAS
# cid.map.cas.pc <- get_cid(
#   clean$cas,
#   from = "cas",
#   match = "all",
#   verbose = FALSE)
#
# saveRDS(cid.map.cas.pc, )
cid.map.cas.pc <- readRDS("cid.map.cas.pc_0608.RDS")

cid.map.cas.pc.uniq <- cid.map.cas.pc %>%
  group_by(query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  rename(cas = query) %>%
  left_join(clean.inch, join_by(cas == cas))

# putting everythign together
# inch.map.cas.cir.uniq -> clean.inch
# cid.map.inch.pc.uniq , cid.map.name.pc.uniq, cid.map.cas.pc.uniq

glimpse(cid.map.inch.pc.uniq)
glimpse(cid.map.name.pc.uniq)
glimpse(cid.map.cas.pc.uniq)

# Combine all cid results, keep if record unique, OR if duplicate records are the same
cid.map <- bind_rows(cid.map.inch.pc.uniq,
                        cid.map.name.pc.uniq,
                        cid.map.cas.pc.uniq) %>%
  group_by(tolower(name)) %>%
  filter(n_distinct(cid) == 1) %>%
  slice(1) %>%
  ungroup()

# VALIDATION
# empty, good
cid.map %>% group_by(name) %>% summarise(n()) %>% filter(`n()`>1)
temp <- cid.map %>% group_by(cid) %>% summarise(n()) %>%
  filter(`n()`>1) %>% ungroup()
# some duplicate cids, many seem to be slight vairations in naming
cid.dupe <- cid.map %>% filter(cid %in% temp$cid) %>% arrange(cid)
# exclude cid.dupe from cid.map
cid.map.uniq <- cid.map %>% filter(!(cid %in% cid.dupe$cid))

# retrieve other data for table Chemicals based on CID
  ids.other <- pc_prop(cid.map.uniq$cid, properties = c("Title", "InChIKey", "InChI",
                                                 "SMILES","IUPACName",
                                                 "MolecularFormula"), verbose = F)
saveRDS(ids.other, "ids.other.RDS")

# retrieve CAS based on CID for verification
cas.pc <- pc_sect(cid.map.uniq$cid, "CAS", domain = "compound", verbose=F)
saveRDS(cas.pc, "cas.pc.RDS")

# Group and take only the uniques
ids.other.uniq <- ids.other %>%
  group_by(CID) %>%
  filter(n() == 1) %>%
  slice(1) %>%
  ungroup()

cas.pc.uniq <- cas.pc %>%
  group_by(CID) %>%
  filter(n() == 1) %>%
  slice(1) %>%
  ungroup()

# combine EVERYTHING, then compare and only keep if values agree (NAs ok)
clean.ids <- clean %>%
  left_join(clean.inch %>% select(name, inchikey), join_by(name == name)) %>%
  left_join(cid.map.uniq %>% select(name, inchikey, cid, cas), join_by(name == name)) %>%
  left_join(cas.pc.uniq %>% select(Name, CID, Result), join_by(cid == CID)) %>%
  left_join(ids.other.uniq, join_by(cid == CID))


# 1. Check that cas.x = cas.y
# if not equal, check cas.x is actually cas
# if not, use cas.y, if yes, then remove entry
# 2. Check inchikey.x = inchikey.y, if not remove entry


clean.pen <- clean.ids %>%
  mutate(
    cas.final = case_when(
      cas.x == cas.y ~ cas.x,
      !is.cas(cas.x) ~ cas.y,
      TRUE ~ NA_character_
    ),
    inchikey.final = apply(
      cbind(inchikey.x, inchikey.y, InChIKey), # bind the 3 values
      1, # eval by row
      function(v) { # v is the current row vector (created by cbind)

        v <- v[!is.na(v)] # remove the na vals

        if (length(v) == 0) return(NA_character_) # if all 3 cols are na, na
        if (length(unique(v)) == 1) return(v[1]) # if 1 unique val, that one!
        NA_character_ # if 2 or 3 uniques, remove :(
      }
    )
  ) %>%
  filter(!is.na(inchikey.final)) %>% # since inchikey is the pk, it cannot be na
  select(name, label, inchikey = inchikey.final, cas = cas.final, cid,
         formula = MolecularFormula, smiles = SMILES, inchi = InChI,
         iupac = IUPACName)

# check how many nas there are for the fields
colSums(is.na(clean.pen))
clean.pen %>% filter(is.na(cid)) %>% select(inchikey)
clean.patch <- clean.pen %>% filter(is.na(cid))

# Get the remaining missing properties...
t <- get_cid(clean.patch$inchikey,
    from = "inchikey",
    match = "all",
    verbose = FALSE)

t2 <- t %>% filter(!is.na(cid)) %>%
  group_by(query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  rename(inchikey = query)

clean.pen.cid <- clean.pen %>% left_join(t2, join_by(inchikey == inchikey)) %>%
  mutate(
    cid = case_when(
      !(is.na(cid.x)) ~ cid.x,
      TRUE ~ cid.y
    )
  ) %>%
  select(-cid.x, -cid.y)


clean.pen.ids <- pc_prop(clean.pen.cid$cid, properties = c("Title", "InChIKey", "InChI",
                                                      "SMILES","IUPACName",
                                                      "MolecularFormula"), verbose = F)

clean.final <- clean.pen.cid %>%
  left_join(clean.pen.ids %>% filter(!is.na(CID)), join_by(cid == CID)) %>%
  group_by(inchikey) %>%
  filter(n() == 1) %>%
  ungroup()

# Both should be true...
nrow(clean.final) == length(unique(clean.final$inchikey))
nrow(clean.final) == length(unique(clean.final$cid)) + sum(is.na(clean.final$cid)) - 1

saveRDS(clean.final, "clean.final.RDS")

# retrieve a descripton for each based on CID
desc <- pc_sect(clean.final$cid, "Record Description", domain = "compound", verbose=F)
saveRDS(desc, "desc.RDS")

glimpse(desc)
view(desc)

# Retrieve identifier info from EPA based on inchikey and name
bpa <- chemical_equal_batch(word_list = clean.all$inchikey,
                            rate_limit = 0.3)

bpa.n <- chemical_equal_batch(word_list = clean$name,
                            rate_limit = 0.3)

# consolidate bpa and bpa.n, keep only when at least 1 id is equal or missing
t <- bpa %>% left_join(clean.all %>% select(inchikey, name),
                       join_by(searchValue == inchikey))

bpa.all <- bpa.n %>% left_join(t, join_by(searchValue == name)) %>%
  filter((dtxcid.x == dtxcid.y | xor(is.na(dtxcid.x), is.na(dtxcid.y))) &
           (dtxsid.x == dtxsid.y | xor(is.na(dtxsid.x), is.na(dtxsid.y))) &
           (smiles.x == smiles.y | xor(is.na(smiles.x), is.na(smiles.y)))) %>%
  mutate(name = searchValue, inchikey = searchValue.y,
         dtxcid = coalesce(dtxcid.x, dtxcid.y),
         dtxsid = coalesce(dtxsid.x, dtxsid.y),
         casrn = coalesce(casrn.x, casrn.y),
         smiles = coalesce(smiles.x, smiles.y),
         epaname = coalesce(preferredName.x, preferredName.y)
         ) %>%
  select(name, inchikey, dtxcid, dtxsid, casrn, smiles, epaname)

# Retrieve chemical details using the InChIKey
epa.details.sid <- get_chemical_details_batch(DTXSID = bpa.all$dtxsid)
epa.details.cid <- get_chemical_details_batch(DTXCID = bpa.all$dtxcid)

t1 <- get_aggregate_records_by_dtxsid(DTXSID = bpa.all$dtxsid[1])
t2 <- get_bioactivity_details(DTXSID = bpa.all$dtxsid[1])
# t3 <- get_bioactivity_summary(DTXSID = bpa.all$dtxsid[1])
t4 <- get_biomonitoring_data(DTXSID = bpa.all$dtxsid[1])
t5 <- get_cancer_hazard(DTXSID = bpa.all$dtxsid[1])
library(png)
t6 <- get_chemical_image(DTXSID = bpa.all$dtxsid[1])
t7 <- get_chem_info(DTXSID = bpa.all$dtxsid[1])
t8 <- get_chem_props_exp(DTXSID = bpa.all$dtxsid[1])
t9 <- get_chem_props_pred(DTXSID = bpa.all$dtxsid[1])
t10 <- get_chem_props_summary(DTXSID = bpa.all$dtxsid[1])
t11 <- get_demographic_exposure_prediction(DTXSID = bpa.all$dtxsid[1])
t12 <- get_exposure_functional_use(DTXSID = bpa.all$dtxsid[1])
t13 <- get_exposure_product_data(DTXSID = bpa.all$dtxsid[1])
t14 <- get_fate_by_dtxsid(DTXSID = bpa.all$dtxsid[1])
t15 <- get_general_exposure_prediction(DTXSID = bpa.all$dtxsid[1])
t16 <- get_general_use_keywords(DTXSID = bpa.all$dtxsid[1])
t17 <- get_genetox_details(DTXSID = bpa.all$dtxsid[1])
t18 <- get_genetox_summary(DTXSID = bpa.all$dtxsid[1])
t19 <- get_hazard_by_dtxsid(DTXSID = bpa.all$dtxsid[1])
t20 <- get_reported_functional_use(DTXSID = bpa.all$dtxsid[1])
t21 <- get_single_sample_records_by_dtxsid(DTXSID = bpa.all$dtxsid[1])
t22 <- get_skin_eye_hazard(DTXSID = bpa.all$dtxsid[1])


##################WAIT FOR API KEY FROM EPA


################do chemspider
# test
get_csid("triclosan")

csids <- get_csid(clean.all$inchikey, from = "inchikey", match = "all")

# cs_convert()
################

desc.comb <- desc %>% group_by(CID) %>%
  summarise(Result= paste(Result, collapse = ". "), .groups = "drop")


clean.all <- clean.final %>%
  left_join(desc.comb, join_by(cid == CID)) %>%
  rename(desc = Result)

write.xlsx(clean.all,
           "CEC_Table_Chemicals_0622.xlsx")





# List that needs manuallly attention
# 1. no unique inchikey or cid based on name or CAS (ie not in clean.ids)
# 2. different inchikey from cas and cid (in inchi.valid.no)
# 3. different cas from provided and cid (in cas.valid.no)

attn <- as_tibble(setdiff(clean$name, clean.all$name)) %>% rename(name = value)
attn.list <- attn %>% left_join(clean, join_by(name == name))

write.xlsx(attn.list %>% rename(`Original Name` = name, `Original CAS` = cas),
           "CEC Manual Verification List.xlsx")

# TABLE: PROPERTIES -------------------------------------------------------

# Physical Properties
color <- pc_sect(clean.all$cid, "Color / Form", domain = "compound", verbose=F)
odor <- pc_sect(clean.all$cid, "Odor", domain = "compound", verbose=F)
bp <- pc_sect(clean.all$cid, "Boiling Point", domain = "compound", verbose=F)
mp <- pc_sect(clean.all$cid, "Melting Point", domain = "compound", verbose=F)
decomposition <- pc_sect(clean.all$cid, "Decomposition", domain = "compound", verbose=F)

saveRDS(pp.colour, "color.RDS")
saveRDS(pp.odor, "odor.RDS")
saveRDS(pp.bp, "bp.RDS")
saveRDS(pp.mp, "mp.RDS")
saveRDS(pp.decompn, "decomposition.RDS")

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

# Chemical Releases

fate.exp <- pc_sect(cid.map.cl2$cid, "Environmental Fate/Exposure Summary",
                domain = "compound", verbose=F)

# Environmental Effects

fate <- pc_sect(cid.map.cl2$cid, "Environmental Fate",
                    domain = "compound", verbose=F)

# FIELD: Human Health Effects, multiple sources
health.eff <- pc_sect(cid.map.cl2$cid, "Health Effects",
                      domain = "compound", verbose=F)
target.org <- pc_sect(cid.map.cl2$cid, "Target Organs",
                      domain = "compound", verbose=F)
adv.eff <- pc_sect(cid.map.cl2$cid, "Adverse Effects",
                      domain = "compound", verbose=F)
signs <- pc_sect(cid.map.cl2$cid, "Signs and Symptoms",
                 domain = "compound", verbose=F)
human.health <- bind_rows(health.eff, target.org, adv.eff, signs) %>% arrange(CID)

# Sources/Uses
sources <- pc_sect(cid.map.cl2$cid, "Uses",
                   domain = "compound", verbose=F) %>% filter(!is.na(Result))

# Retain 1st result of SourceName for "Chemical Use (Source)"
# "Haz-Map, Information on Hazardous Chemicals and Occupational Diseases"
# Chemical Use (souce)
chem.use <- sources %>%
  filter(SourceName=="Haz-Map, Information on Hazardous Chemicals and Occupational Diseases") %>%
  group_by(CID) %>% summarise(Name = first(Name),
                              Result = first(Result))



# Source Industry
industry <- sources %>%
  group_by(CID) %>% summarise(Name = first(Name),
                              Result =
                                Result[grepl("Category: Industry", Result,
                                             ignore.case = TRUE)][1])





# SOURCE: Effluent Concentrations is likely more suitable
eff.conc <- pc_sect(cid.map.cl2$cid, "Effluent Concentrations",
                domain = "compound", verbose=F) %>% filter(!is.na(Result))

# SOURCE: Environmental Water Concentrations
env.water <- pc_sect(cid.map.cl2$cid, "Environmental Water Concentrations",
                 domain = "compound", verbose=F) %>% filter(!is.na(Result))

# FIELD: Drinking Water Sources
drink <- water %>%
  filter(str_detect(Result, regex(paste0(
    "(drinking water|tap water|finished water|treated water|",
    "source water|raw water|intake|distribution system|",
    "water utility|potable water)"), ignore_case = TRUE)))

# Sources Facility/ Location
facility <- water %>%
  filter(!(Result %in% drink$Result), str_detect(Result,regex(paste0(
    #not sure if ^ should be included to force mutual exclusivity
    # --- Direct facility + discharge ---
    "(",
    "(plant|facility|factory|industry|industrial|manufacturer|mill|refinery|site|",
    "landfill|leachate|wastewater treatment plant|wwtp|sewage treatment plant)",
    ".*",
    "(effluent|influent|discharge|release|outfall)",
    ")",
    "|",
    # --- Indirect attribution (spatial linkage) ---
    "((downstream|plume|near)\\s+(of\\s+)?",
    "(plant|facility|factory|industry|landfill|wwtp|sewage treatment plant)",
    ")"
  ),
  ignore_case = TRUE)))


# FIELD: Receiving Watersheds
watershed <- water %>%
  filter(!(Result %in% drink$Result), !(Result %in% facility$Result),
         #not sure if ^ should be included to force mutual exclusivity
    str_detect(
      Result,
      regex(
        paste0(
          "(",
          # core environmental media
          "river|stream|lake|surface water|groundwater|aquifer|basin|watershed",
          "|",
          # atmospheric deposition media (C-specific case)
          "rain|rainwater|snow|fog",
          ")",
          ".*",
          "(",
          # distribution / monitoring signals
          "across|multiple|various|survey|monitoring|samples|sites|regions|areas",
          "|",
          # quantitative monitoring language (strong signal)
          "detected in|frequency of detection|%|concentration|ng/L|ug/L",
          ")"
        ),
        ignore_case = TRUE
      )
    )
  ) %>%
  filter(
    !str_detect(
      Result,
      regex(
        paste0(
          # exclude source attribution (B)
          "plant|facility|wwtp|wastewater|landfill|effluent|discharge|outfall",
          "|",
          # exclude pathway / mechanism (A)
          "runoff|enter|transport|leach|deposition into|washoff|scavenging|downstream|plume"
        ),
        ignore_case = TRUE
      )
    )
  )


# FIELD: Chemical Entry Points
# SOUCRE: Environmental Fate / Exposure Summary
entry <- pc_sect(cid.map.cl2$cid, "Environmental Fate / Exposure Summary",
               domain = "compound", verbose=F)

# Monitoring Requirements
req <- pc_sect(cid.map.cl2$cid, "Regulatory Information",
                 domain = "compound", verbose=F)

# Removal Technologies
tech <- pc_sect(cid.map.cl2$cid, "Environmental Biodegradation",
               domain = "compound", verbose=F)

# Safe Production
# safe <- pc_sect(cid.map.cl2$cid, "Storage Conditions",
#                 domain = "compound", verbose=F)
# Safe Use
uses <- pc_sect(cid.map.cl2$cid, "Preventive Measures",
               domain = "compound", verbose=F)

# Safe Disposal
disposal <- pc_sect(cid.map.cl2$cid, "Disposal Methods",
               domain = "compound", verbose=F)

# Consumer Products
# consum <- pc_sect(head(cid.map.cl2)$cid, "Consumer Uses",
#                     domain = "compound", verbose=F)

# Exposure Routes

exp.routes <-  pc_sect(cid.map.cl2$cid, "Exposure Routes",
                          domain = "compound", verbose=F)

# Exposure Baseline

# exp.base <- pc_sect(cid.map.cl2$cid, "Metabolism / Metabolites",
#                           domain = "compound", verbose=F)
# 8.3 Metabolism / Metabolites (P2)


# Transgenerational Effects

# generation <- pc_sect(cid.map.cl2$cid, "Health Effects",
#                         domain = "compound", verbose=F)

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
openxlsx::write.xlsx(merged.ordered, "cec_pipeline_data_review_v1.xlsx")
# writexl::write_xlsx(merged.pc.name, "write.xlsx")



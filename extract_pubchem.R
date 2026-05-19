library(readxl)
library(tidyverse)
library(webchem)
library(janitor)

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

# CIR API -----------------------------------------------------------------
# Get InChIKey from CAS, CIR api
inch.map.cas.cir <- cir_query(clean$CAS,
                               representation = "stdinchikey",
                               match = "all",
                               verbose = FALSE)

inch.map.cas.cir.uniq <- inch.map.cas.cir %>%
  group_by(CAS = query) %>%
  filter(n()==1, !(is.na(stdinchikey))) %>% ungroup() %>%
  left_join(clean %>%
              select(`Chemical Name`, CAS),
            join_by(query == CAS)) %>%
  rename(InChIKey = stdinchikey) %>%
  reframe(`Chemical Name`, CAS, InChIKey) %>%
  mutate(InChIKey = str_replace(InChIKey, "InChIKey=", ""))


# PC API ------------------------------------------------------------------
# Get CID from inchikey (retrieved through CIR)
cid.map.inch.pc <- get_cid(
  inch.map.cas.cir.uniq$InChIKey,
  from = "inchikey",
  match = "all",
  verbose = FALSE)

cid.map.inch.pc.uniq <- cid.map.inch.pc %>%
  group_by(InChIKey = query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  left_join(inch.map.cas.cir.uniq,
            join_by(InChIKey == InChIKey)) %>%
  rename(CID = cid) %>%
  reframe(`Chemical Name`, CAS, CID)


# Get CID from chemical name
cid.map.name.pc <- get_cid(
  clean$`Chemical Name`,
  from = "name",
  match = "all",
  verbose = FALSE)

cid.map.name.pc.uniq <- cid.map.name.pc %>%
  group_by(`Chemical Name` = query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  left_join(clean %>%
              select(`Chemical Name`, CAS),
            join_by(query == `Chemical Name`)) %>%
  rename(CID = cid) %>%
  reframe(`Chemical Name`, CAS, CID)

# Get CID from CAS
cid.map.cas.pc <- get_cid(
  clean$CAS,
  from = "cas",
  match = "all",
  verbose = FALSE)

cid.map.cas.pc.uniq <- cid.map.cas.pc %>%
  group_by(query) %>%
  filter(n()==1, !(is.na(cid))) %>% ungroup() %>%
  left_join(clean %>%
              select(`Chemical Name`, CAS),
            join_by(query == CAS)) %>%
  rename(CID = cid, CAS = query) %>%
  reframe(`Chemical Name`, CAS, CID)

# keep if record unique, OR if duplicate records are the same
cid.map <- bind_rows(cid.map.inch.pc.uniq,
                        cid.map.name.pc.uniq,
                        cid.map.cas.pc.uniq) %>%
  group_by(`Chemical Name`) %>%
  filter(n() == 1 | (n() > 1 & n_distinct(across(everything())) == 1)) %>%
  slice(1) %>% ungroup() %>%
  select(`Chemical Name`, CID)
#
# cid.map.dupes <- get_dupes(cid.map, CID)[,1] %>% mutate(CID = as.numeric(CID))
cid.map.uniq <- cid.map %>%
  group_by(CID) %>%
  filter(n() == 1) %>%
  ungroup()

attn <- clean %>%
  filter(!(`Chemical Name` %in% cid.map.uniq$`Chemical Name`)) %>%
  select(`Chemical Name`)

# Retrieve Column Data ----------------------------------------------------

# Add CAS to cid.map
cas.db <- clean %>% filter(is.cas(clean$CAS)) %>% select(`Chemical Name`, CAS)
cas.pc <- pc_sect(cid.map.uniq$CID, section = "CAS", verbose = FALSE) %>% filter(!is.na(Result))

# 3 = name, cid, cas
cid.map3 <- cid.map.uniq %>%
  left_join(cas.pc, join_by(CID == CID))

# Add InChIKey to cid.map



#
id.all <- pc_prop(cid.map$CID, properties = c("InChIKey","InChI",
                                                       "SMILES","IUPACName",
                                                       "MolecularFormula"
                                                       ), verbose = F) %>%
  left_join(cid.map %>% select(`Chemical Name`, CID), join_by(CID==CID))

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



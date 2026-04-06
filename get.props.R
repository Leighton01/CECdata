library(readxl)
library(tidyverse)
library(webchem)

# Ping to see if api works
ping_service(
  service = c("bcpc", "chebi", "chembl", "cs", "cs_web", "cir", "cts", "etox", "fn",
              "nist", "opsin", "pc", "srs", "wd"),
  apikey = NULL
)

load(cid.map.cl2)


# Function: Group and Save RDS --------------------------------------------

group.save <- function(df){

  df.comb <- df %>% group_by(CID) %>%
    summarise(Result = paste(Result, collapse = ".\n"), .groups = "drop")

  df.comb <- df.comb %>% left_join(cid.map.cl, join_by (CID == cid))

  saveRDS(df.comb, paste(get_name(df), ".RDS"))

  write.csv(df.comb, paste(get_name(df), ".csv"),
            fileEncoding="Windows-1252", row.names = FALSE)

}



# Chemical Properties -----------------------------------------------------

desc1 <- pc_sect(cid.map.cl2$cid, "Agrochemical Information", domain = "compound", verbose=F)
desc2 <- pc_sect(cid.map.cl2$cid, "Biologic Information", domain = "compound", verbose=F)

# Physical Properties -----------------------------------------------------
pp.colour <- pc_sect(cid.map.cl2$cid, "Color / Form", domain = "compound", verbose=F)
pp.odor <- pc_sect(cid.map.cl2$cid, "Odor", domain = "compound", verbose=F)
pp.bp <- pc_sect(cid.map.cl2$cid, "Boiling Point", domain = "compound", verbose=F)
pp.mp <- pc_sect(cid.map.cl2$cid, "Melting Point", domain = "compound", verbose=F)
pp.decompn <- pc_sect(cid.map.cl2$cid, "Decomposition", domain = "compound", verbose=F)

saveRDS(pp.colour, file="pp.colour.RDS")
saveRDS(pp.odor, file="pp.odor.RDS")
saveRDS(pp.bp, file="pp.bp.RDS")
saveRDS(pp.mp, file="pp.mp.RDS")
saveRDS(pp.decompn, file="pp.decompn.RDS")

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

saveRDS(pp.comb, "pp.comb.RDS")
#
# pp.by.name <- pp.comb %>% left_join(cid.map.cl, join_by (CID == cid)) %>% glimpse
# write.csv(pp.by.name, "pp.name.csv", fileEncoding="Windows-1252", row.names = FALSE)



# Chemical Releases -------------------------------------------------------
fate.exp <- pc_sect(cid.map.cl2$cid, "Environmental Fate/Exposure Summary",
                domain = "compound", verbose=F)
# One result per CID

saveRDS(fate.exp, "fate.exp.RDS")


# Environmental Effects ---------------------------------------------------
fate <- pc_sect(cid.map.cl2$cid, "Environmental Fate",
                    domain = "compound", verbose=F)

saveRDS(fate, "fate.RDS")

# Human Health Effects ----------------------------------------------------
signs <- pc_sect(cid.map.cl2$cid, "Signs and Symptoms",
                 domain = "compound", verbose=F)
saveRDS(signs, "signs.RDS")

# Sources -----------------------------------------------------------------
# Chemical Use (Source)
# Sources Facility / Location
# Source Industry







# Water Sources and Watersheds ------------------------------------------
water <- pc_sect(cid.map.cl2$cid, "Environmental Water Concentrations",
        domain = "compound", verbose=F)
saveRDS(water, "water.RDS")

# Drinking Water Sources
drink <-

# Receiving Watersheds, paragraph 2
shed <-


# Chemical Entry Points ---------------------------------------------------

12.2.8 Environmental Fate (S2P1)


# Monitoring Requirements -------------------------------------------------
req <- pc_sect(cid.map.cl2$cid, "Regulatory Information",
                 domain = "compound", verbose=F)

saveRDS(req, "req.RDS")
# Removal Technologies ----------------------------------------------------
tech <- pc_sect(cid.map.cl2$cid, "Environmental Biodegradation",
               domain = "compound", verbose=F)

saveRDS(tech, "tech.RDS")
# Safe Production ---------------------------------------------------------
safe <- pc_sect(cid.map.cl2$cid, "Storage Conditions",
                domain = "compound", verbose=F)
saveRDS(safe, "safe.RDS")
# Safe Use ----------------------------------------------------------------
use <- pc_sect(cid.map.cl2$cid, "Personal Protective Equipment",
               domain = "compound", verbose=F)
saveRDS(use, "use.RDS")
11.5.2 Personal Protective Equipment (PPE)


# Safe Disposal -----------------------------------------------------------


11.3.2 Disposal Methods (P4(

# Consumer Products -------------------------------------------------------


  9.1.1 Household Products




# Exposure Routes ---------------------------------------------------------


12.1.9 Exposure Routes


# Exposure Baseline -------------------------------------------------------


8.3 Metabolism / Metabolites (P2)


# Transgenerational Effects -----------------------------------------------



12.1.8 Health Effects

# Hormetic Effects --------------------------------------------------------


12.1.6 Evidence for Carcinogenicity



# H2O Sol. ----------------------------------------------------------------


12.2.13 Volatilization from Water / Soil
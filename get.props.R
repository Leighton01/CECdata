library(readxl)
library(tidyverse)
library(webchem)

# Ping to see if api works
ping_service(
  service = c("bcpc", "chebi", "chembl", "cs", "cs_web", "cir", "cts", "etox", "fn",
              "nist", "opsin", "pc", "srs", "wd"),
  apikey = NULL
)

# load(cid.map.cl2)


# # Function: Group and Save RDS (Deprecated) --------------------------------------------
# # : compiling a package instead

# Retrieve Column Data ----------------------------------------------------
## Chemical Properties ----------------------------------------------------

# desc1 <- pc_sect(cid.map.cl2$cid, "Agrochemical Information", domain = "compound", verbose=F)
desc2 <- pc_sect(cid.map.cl2$cid, "Record Description", domain = "compound", verbose=F)

# saveRDS(desc1, "desc1.RDS")
# saveRDS(desc2, "desc2.RDS")


## Physical Properties -----------------------------------------------------
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

# saveRDS(pp.comb, "pp.comb.RDS")
#
# pp.by.name <- pp.comb %>% left_join(cid.map.cl, join_by (CID == cid)) %>% glimpse
# write.csv(pp.by.name, "pp.name.csv", fileEncoding="Windows-1252", row.names = FALSE)



## Chemical Releases -------------------------------------------------------
fate.exp <- pc_sect(cid.map.cl2$cid, "Environmental Fate/Exposure Summary",
                domain = "compound", verbose=F)


## Environmental Effects --------------------------------------------------
fate <- pc_sect(cid.map.cl2$cid, "Environmental Fate",
                    domain = "compound", verbose=F)

# group.save(fate, "fate.RDS")

## Human Health Effects ----------------------------------------------------
signs <- pc_sect(cid.map.cl2$cid, "Signs and Symptoms",
                 domain = "compound", verbose=F)
# group.save(signs, "signs.RDS")

## Sources -----------------------------------------------------------------
sources <- pc_sect(cid.map.cl2$cid, "Sources/Uses",
                   domain = "compound", verbose=F)

# Chemical Use (Source)
# Sources Facility / Location
# Source Industry

## Water Sources and Watersheds ------------------------------------------
water <- pc_sect(cid.map.cl2$cid, "Environmental Water Concentrations",
        domain = "compound", verbose=F)
# group.save(water)

# Drinking Water Sources
# drink <-

# Receiving Watersheds, paragraph 2
# shed <-


## Chemical Entry Points ---------------------------------------------------

# 12.2.8 Environmental Fate (S2P1)


## Monitoring Requirements -------------------------------------------------
req <- pc_sect(cid.map.cl2$cid, "Regulatory Information",
                 domain = "compound", verbose=F)

## Removal Technologies ----------------------------------------------------
tech <- pc_sect(cid.map.cl2$cid, "Environmental Biodegradation",
               domain = "compound", verbose=F)

## Safe Production ---------------------------------------------------------
safe <- pc_sect(cid.map.cl2$cid, "Storage Conditions",
                domain = "compound", verbose=F)
## Safe Use ----------------------------------------------------------------
uses <- pc_sect(cid.map.cl2$cid, "Personal Protective Equipment",
               domain = "compound", verbose=F)

## Safe Disposal -----------------------------------------------------------
disposal <- pc_sect(cid.map.cl2$cid, "Disposal Methods",
               domain = "compound", verbose=F)

## Consumer Products -------------------------------------------------------

consum <- pc_sect(cid.map.cl2$cid, "Household Products",
                    domain = "compound", verbose=F)


## Exposure Routes ---------------------------------------------------------

exp.routes <-  pc_sect(cid.map.cl2$cid, "Exposure Routes",
                          domain = "compound", verbose=F)
# 12.1.9 Exposure Routes


## Exposure Baseline -------------------------------------------------------

exp.base <- pc_sect(cid.map.cl2$cid, "Metabolism / Metabolites",
                          domain = "compound", verbose=F)
# 8.3 Metabolism / Metabolites (P2)


## Transgenerational Effects -----------------------------------------------

generation <- pc_sect(cid.map.cl2$cid, "Health Effects",
                        domain = "compound", verbose=F)

# 12.1.8 Health Effects

## Hormetic Effects --------------------------------------------------------

carcin <- pc_sect(cid.map.cl2$cid, "Evidence for Carcinogenicity",
                      domain = "compound", verbose=F)

## H2O Sol. ----------------------------------------------------------------

h20 <- pc_sect(cid.map.cl2$cid, "Volatilization from Water / Soil",
                      domain = "compound", verbose=F)



# Combine -----------------------------------------------------------------
list1 <- lst( desc2,
              pp.comb,
              fate.exp,
              fate,
              signs,
              water,
              req,
              tech,
              safe,
              # uses, # no result
              disposal,
              consum,
              exp.routes,
              exp.base,
              generation,
              carcin,
              h20)

# check if any is empty (prob shouldve checked earlier)
which((vapply(list1, function(x) !("Result" %in% colnames(x)), logical(1))))

# Drop unneeded columns
list2 <- lst()

# Clean
for (i in seq_along(list1)){

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

# Add name
merged.pc.name <- merged.pc.df %>% left_join(cid.map.cl, join_by (CID == cid))


# write.csv(merged.pc.name, "merged.pc.df.csv",
#           fileEncoding="UTF-8", row.names = FALSE)

openxlsx::write.xlsx(merged.pc.name, "open.xlsx")
# writexl::write_xlsx(merged.pc.name, "write.xlsx")



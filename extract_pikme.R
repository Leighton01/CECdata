# Downloaded PikMe from: https://zenodo.org/records/15647470
# Downloaded on 20260511


# unforutantely does not contain the info we need, might use vannmiljo_paramater_id_norway to get MEC from vannmiljo db but will try the norman ones first
library(tidyverse)
library(arrow)
library(httr2)
library(jsonlite)
library(Vannmiljo-R)

raw.pikme <- read_parquet("data/pikme_all.parquet")
# count overlap, pretty good!
sum(unique(raw.pikme$inchikey) %in% t.chemicals$inchikey)

van.ids <- raw.pikme %>% filter(!is.na(vannmiljo_paramater_id_norway),inchikey %in% t.chemicals$inchikey) %>%
  select(vannmiljo_paramater_id_norway) %>% distinct()



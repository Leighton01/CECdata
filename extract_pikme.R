# Downloaded PikMe from: https://zenodo.org/records/15647470
# Downloaded on 20260511

library(tidyverse)
library(arrow)

raw.pikme <- read_parquet("data/pikme_all.parquet")
# count overlap, pretty good!
sum(unique(raw.pikme$inchikey) %in% merged.ordered$InChIKey)


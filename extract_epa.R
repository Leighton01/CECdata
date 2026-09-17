library(readxl)
library(openxlsx)
library(tidyverse)
library(webchem)
library(janitor)
library(ctxR)

# register epa api key, store in envir
register_ctx_api_key("4d4cbbce-4265-4b24-9712-326c0d8ea618", write = T)

# get ids
t.chemicals <- readRDS("t.chemicals.RDS")
epa.info <- readRDS("epa.info.RDS")

# get what the chemicals were used for
# func.use <- get_exposure_functional_use_batch(DTXSID=unique(t.chemicals$dtxsid))
# prod.data <- get_exposure_product_data_batch(DTXSID=unique(t.chemicals$dtxsid))
# saveRDS(func.use, "func.use.RDS")
# saveRDS(prod.data, "prod.data.RDS")
readRDS("prod.data.RDS")

# we just need prod.data
prod.data[[1]] %>% colnames()
prod.data[[1]]$prodfam %>% unique()
prod.data[[1]]$prodtype %>% unique()
prod.data[[1]]$productname %>% unique()


prod.data.df <- dplyr::bind_rows(lapply(prod.data, function(x) {
  data.frame(
    dtxsid = x$dtxsid,
    prodfam = x$prodfam,
    prodtype = x$prodtype,
    productname = x$productname
  )
})) %>%
  distinct() %>%
  filter(!is.na(dtxsid))

# add inchikey to it,
temp <- prod.data.df %>% left_join(t.chemicals
                                   %>% select(dtxsid, inchikey),
                                   join_by(dtxsid == dtxsid)) %>%
  filter(!is.na(inchikey)) %>%
  distinct() %>%
  select(inchikey, prodfam, prodtype, productname)

prod.data.long <- temp %>%
  pivot_longer(cols = c(prodfam, prodtype, productname),
               names_to = "info_type",
               values_to = "content") %>%
  mutate(pk = paste0(inchikey, "_", info_type))

prod.data.long %>% glimpse()
t.info3 <- t.info2 %>% select(-CID) %>% rbind(prod.data.long)
openxlsx::write.xlsx(t.info3, "CEC_Table_Information_20260917.xlsx")

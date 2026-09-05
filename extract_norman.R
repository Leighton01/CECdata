library(httr2)
library(dplyr)
library(purrr)
library(tidyr)

# get info from norman suspect exchange empodat db

# write function to get data by inchikey
get_empodat <- function(inchikey) {
  purrr::map(
    inchikey,
    \(key) {
      request(
        paste0(
          "https://www.norman-network.com/nds/api/empodat/inchikey/",
          URLencode(key, reserved = TRUE),
          "/JSON"
        )
      ) %>%
        req_perform() %>%
        resp_body_json(simplifyVector = TRUE)
    }
  )
}

# extract

empodat <- get_empodat(unique(t.chemicals$inchikey))

# empodat <- t.chemicals %>%
#   distinct(inchikey) %>%
#   filter(!is.na(inchikey), inchikey != "") %>%
#   mutate(data = map(inchikey, get_empodat)) %>%
#   unnest(data)

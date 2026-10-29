library(terra); library(dplyr)

files   <- list.files("data_cache/lst", full.names = TRUE)
got     <- as.Date(sub(".*lst_(.*)\\.tif", "\\1", files))
n_valid <- sapply(files, \(f) global(!is.na(rast(f)), "sum")[1, 1])

# MODIS 8-day blocks always start on 1 January, then every 8 days: 46 per year
expected <- do.call(c, lapply(2021:2023, \(y)
                              seq(as.Date(paste0(y, "-01-01")), by = 8, length.out = 46)))

inv <- data.frame(date = expected) |>
  left_join(data.frame(date = got, n_valid = unname(n_valid)), by = "date") |>
  mutate(n_valid   = coalesce(n_valid, 0),
         pct_valid = round(100 * n_valid / 650, 1),
         status = case_when(!date %in% got ~ "missing composite",
                            n_valid == 0   ~ "empty (all cloud)",
                            pct_valid < 50 ~ "mostly cloud",
                            TRUE           ~ "usable"))

table(inv$status)
write.csv(inv, "data_cache/lst_inventory.csv", row.names = FALSE)

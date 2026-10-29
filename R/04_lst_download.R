library(rstac); library(terra); library(sf)

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
bbox <- as.numeric(st_bbox(st_transform(city, 4326)))
pc   <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
dir.create("data_cache/lst", showWarnings = FALSE)

# Quality rule: keep good pixels (bits 0-1 = 0),
# or "other quality" pixels (bits 0-1 = 1) with error ≤ 2 K (bits 6-7 ≤ 1)
qc_keep <- function(qc) {
  mand <- qc %% 4
  err  <- (qc %/% 64) %% 4
  (mand == 0) | (mand == 1 & err <= 1)
}

log <- list()
for (yr in 2021:2023) {
  it <- pc |>
    stac_search(collections = "modis-11A2-061", bbox = bbox,
                datetime = paste0(yr, "-01-01/", yr, "-12-31"), limit = 250) |>
    post_request() |> items_fetch() |>
    items_sign(sign_fn = sign_planetary_computer())
  
  # Terra only (IDs start with MOD)
  it$features <- it$features[startsWith(sapply(it$features, \(f) f$id), "MOD")]
  dates <- as.Date(substr(sapply(it$features, \(f) f$properties$start_datetime), 1, 10))
  message(yr, ": ", length(dates), " Terra composites found")
  
  for (d in as.character(sort(unique(dates)))) {
    out <- sprintf("data_cache/lst/lst_%s.tif", d)
    if (file.exists(out)) next                      # already done, skip it
    
    i <- which(dates == as.Date(d))[1]
    ok <- try({
      area <- vect(st_transform(st_buffer(city, 2000), "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"))
      lst <- crop(rast(paste0("/vsicurl/", it$features[[i]]$assets$LST_Day_1km$href)), area)
      qc  <- crop(rast(paste0("/vsicurl/", it$features[[i]]$assets$QC_Day$href)), area)
      
      lst[lst == 0] <- NA                                   # 0 = no data
      if (isTRUE(global(lst, "max", na.rm = TRUE)[1, 1] > 1000)) lst <- lst * 0.02
      lst <- lst - 273.15                                   # Kelvin → °C
      lst[!qc_keep(qc)] <- NA                               # drop unreliable pixels
      names(lst) <- d
      writeRaster(lst, out, overwrite = TRUE)
      
      n_valid <- global(!is.na(lst), "sum")[1, 1]
      log[[d]] <- data.frame(date = d, n_valid = n_valid, n_cells = ncell(lst))
      message(d, "  valid pixels: ", n_valid, " / ", ncell(lst))
    }, silent = TRUE)
    if (inherits(ok, "try-error")) message(d, "  FAILED, rerun the script later to retry")
  }
}

if (length(log)) {
  new <- do.call(rbind, log)
  old <- if (file.exists("data_cache/lst_log.csv")) read.csv("data_cache/lst_log.csv") else NULL
  write.csv(unique(rbind(old, new)), "data_cache/lst_log.csv", row.names = FALSE)
}
message("DONE. Files on disk: ", length(list.files("data_cache/lst")))

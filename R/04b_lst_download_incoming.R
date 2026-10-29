# ============================================================
# 04b_lst_download_incoming.R - v3: new MODIS weeks (2024) into staging
# v3: full RFC3339 timestamps (Planetary Computer API change)
# ============================================================
library(rstac); library(terra); library(sf)

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
bbox <- as.numeric(st_bbox(st_transform(city, 4326)))
pc   <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
dir.create("data_cache/lst_incoming", showWarnings = FALSE)
sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
area <- vect(st_transform(st_buffer(city, 2000), sinu))      # same area as before

# Same quality rule as the original download (must stay identical)
qc_keep <- function(qc) {
  mand <- qc %% 4
  err  <- (qc %/% 64) %% 4
  (mand == 0) | (mand == 1 & err <= 1)
}

# Full RFC3339 time range, e.g. "2024-01-01T00:00:00Z/2024-01-31T23:59:59Z"
rfc <- function(start, end) paste0(start, "T00:00:00Z/", end, "T23:59:59Z")

# Try a step up to 4 times; return NULL if it keeps failing
retry <- function(f, what, tries = 4) {
  for (i in 1:tries) {
    r <- try(f(), silent = TRUE)
    if (!inherits(r, "try-error")) return(r)
    message(what, ": try ", i, " failed -> ", trimws(as.character(r)))
    Sys.sleep(20)
  }
  NULL
}

log <- list()
for (m in 1:12) {
  start <- as.Date(sprintf("2024-%02d-01", m))
  end   <- seq(start, by = "month", length.out = 2)[2] - 1
  lab   <- format(start, "%Y-%m")
  
  # 1. Search one month only
  it <- retry(function() pc |>
                stac_search(collections = "modis-11A2-061", bbox = bbox,
                            datetime = rfc(start, end), limit = 100) |>
                post_request() |>
                items_sign(sign_fn = sign_planetary_computer()),
              paste(lab, "search"))
  if (is.null(it)) { message(lab, ": SEARCH FAILED, rerun the script later"); next }
  
  # 2. Terra only, and only composites that START in this month (no double counting)
  it$features <- it$features[startsWith(sapply(it$features, \(f) f$id), "MOD")]
  if (length(it$features) == 0) { message(lab, ": no Terra composites"); next }
  dates <- as.Date(substr(sapply(it$features, \(f) f$properties$start_datetime), 1, 10))
  keep  <- dates >= start & dates <= end
  it$features <- it$features[keep]; dates <- dates[keep]
  
  # 3. Download each composite
  for (d in as.character(sort(unique(dates)))) {
    out <- sprintf("data_cache/lst_incoming/lst_%s.tif", d)
    if (file.exists(out) || file.exists(sprintf("data_cache/lst/lst_%s.tif", d))) next
    i <- which(dates == as.Date(d))[1]
    ok <- retry(function() {
      lst <- crop(rast(paste0("/vsicurl/", it$features[[i]]$assets$LST_Day_1km$href)), area)
      qc  <- crop(rast(paste0("/vsicurl/", it$features[[i]]$assets$QC_Day$href)), area)
      lst[lst == 0] <- NA
      if (isTRUE(global(lst, "max", na.rm = TRUE)[1, 1] > 1000)) lst <- lst * 0.02
      lst <- lst - 273.15
      lst[!qc_keep(qc)] <- NA
      names(lst) <- d
      writeRaster(lst, out, overwrite = TRUE)
      lst
    }, paste(d, "download"))
    if (is.null(ok)) { message(d, "  FAILED, rerun later to retry"); next }
    n_valid <- global(!is.na(ok), "sum")[1, 1]
    log[[d]] <- data.frame(date = d, n_valid = n_valid, n_cells = ncell(ok))
    message(d, "  valid pixels: ", n_valid, " / ", ncell(ok))
  }
}

# 4. Log (adds to any earlier log, so reruns don't lose anything)
if (length(log)) {
  new <- do.call(rbind, log)
  old <- if (file.exists("data_cache/lst_incoming_log.csv")) read.csv("data_cache/lst_incoming_log.csv") else NULL
  write.csv(unique(rbind(old, new)), "data_cache/lst_incoming_log.csv", row.names = FALSE)
}
message("DONE. Files in staging: ", length(list.files("data_cache/lst_incoming")))

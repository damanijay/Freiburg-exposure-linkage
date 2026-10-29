# ============================================================
# 10_station_validation.R - satellite LST vs DWD Freiburg station
# ============================================================
library(sf); library(dplyr); library(ggplot2)
dir.create("data_cache/dwd", showWarnings = FALSE)

# ---------- Section 1: download station 01443 (Freiburg), daily climate ----------
base    <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/"
listing <- readLines(base, warn = FALSE)
fname   <- regmatches(listing, regexpr("tageswerte_KL_01443_[0-9_]+_hist\\.zip", listing))[1]
cat("File:", fname, "\n")                           # e.g. tageswerte_KL_01443_19490101_20251231_hist.zip
zipf <- file.path("data_cache/dwd", fname)
if (!file.exists(zipf)) download.file(paste0(base, fname), zipf, mode = "wb")
unzip(zipf, exdir = "data_cache/dwd", overwrite = TRUE)

# ---------- Section 2: station location (from DWD's own metadata) ----------
geo_file <- list.files("data_cache/dwd", pattern = "^Metadaten_Geographie", full.names = TRUE)[1]
geo <- read.table(geo_file, sep = ";", header = TRUE, strip.white = TRUE,
                  fileEncoding = "latin1", quote = "")
loc <- tail(geo, 1)                                 # most recent position (stations move)
stn_elev <- as.numeric(loc[[2]]); lat <- as.numeric(loc[[3]]); lon <- as.numeric(loc[[4]])
cat("Station:", loc[[7]], " lat", lat, " lon", lon, " elevation", stn_elev, "m\n")

hex <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
stn <- st_sf(geometry = st_sfc(st_point(c(lon, lat)), crs = 4326)) |> st_transform(25832)
stn_hex <- st_join(stn, hex["hex_id"])$hex_id
cat("Station lies in hexagon:", stn_hex, "\n")
stopifnot(!is.na(stn_hex))

# ---------- Section 3: daily air temperature 2021-2023 ----------
prod <- list.files("data_cache/dwd", pattern = "^produkt_klima_tag", full.names = TRUE)[1]
dwd  <- read.table(prod, sep = ";", header = TRUE, strip.white = TRUE)
names(dwd) <- trimws(names(dwd))
dwd <- dwd |>
  transmute(day   = as.Date(as.character(MESS_DATUM), "%Y%m%d"),
            tmax  = ifelse(TXK <= -999, NA, TXK),       # daily maximum air temperature
            tmean = ifelse(TMK <= -999, NA, TMK)) |>     # daily mean air temperature
  filter(day >= as.Date("2021-01-01"), day <= as.Date("2023-12-31"))
cat("Station days:", nrow(dwd), " missing tmax:", sum(is.na(dwd$tmax)), "\n")   # expect 1095, ~0

# ---------- Section 4: aggregate station to the same 8-day blocks ----------
inv <- read.csv("data_cache/lst_inventory.csv") |> mutate(date = as.Date(date))
air8 <- inv |>
  transmute(date, end = pmin(date + 7, as.Date(paste0(format(date, "%Y"), "-12-31")))) |>
  rowwise() |>
  mutate(tmax  = mean(dwd$tmax [dwd$day >= date & dwd$day <= end], na.rm = TRUE),
         tmean = mean(dwd$tmean[dwd$day >= date & dwd$day <= end], na.rm = TRUE)) |>
  ungroup()

filled <- read.csv("data_cache/lst_hex_filled.csv") |> mutate(date = as.Date(date))
v <- filled |> filter(hex_id == stn_hex) |> left_join(air8, by = "date") |>
  mutate(month = format(date, "%m"),
         season = ifelse(as.integer(month) %in% c(10:12, 1:3), "cold (Oct-Mar)", "warm (Apr-Sep)"))
print(table(v$source))                              # observed vs filled weeks at the station

# ---------- Section 5: does the satellite track real weather? ----------
obs <- filter(v, source == "observed")
cat("Raw correlation LST vs air Tmax (observed weeks):", round(cor(obs$lst_filled, obs$tmax), 2), "\n")

# Remove the seasonal cycle: anomaly = value minus that month's average
v <- v |> group_by(month) |>
  mutate(lst_anom  = lst_filled - mean(lst_filled[source == "observed"]),
         tmax_anom = tmax - mean(tmax, na.rm = TRUE)) |> ungroup()
obs <- filter(v, source == "observed")
cat("Anomaly correlation (seasons removed):", round(cor(obs$lst_anom, obs$tmax_anom), 2), "\n")

v |> filter(source == "observed") |> group_by(season) |>
  summarise(mean_LST_minus_Tmax = round(mean(lst_filled - tmax), 1), n_weeks = n()) |> print()

# ---------- Section 6: are the FILLED weeks realistic? ----------
fit <- lm(lst_filled ~ tmax, data = obs)            # relationship learned on observed weeks only
v$resid <- v$lst_filled - predict(fit, newdata = v) # how far each week sits from that line
v |> group_by(source) |>
  summarise(n = n(), mean_resid = round(mean(resid), 2),
            rmse_resid = round(sqrt(mean(resid^2)), 2)) |> print()

# ---------- Section 7: figures ----------
p10 <- ggplot(v, aes(tmax, lst_filled, colour = source)) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], linetype = 2) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = paste("Satellite LST vs DWD air temperature, hexagon", stn_hex),
       subtitle = "Dashed line fitted on observed weeks only",
       x = "Station air Tmax, 8-day mean (°C)", y = "MODIS LST (°C)") + theme_minimal()
print(p10); ggsave("outputs/figures/F10_station_scatter.png", p10, width = 7, height = 5, dpi = 200)

p11 <- ggplot(v, aes(date)) +
  geom_line(aes(y = tmax, linetype = "Station air Tmax")) +
  geom_point(aes(y = lst_filled, colour = source), size = 1.6) +
  labs(title = "MODIS LST and DWD air temperature, Freiburg 2021-2023",
       y = "°C", x = NULL, linetype = NULL) + theme_minimal()
print(p11); ggsave("outputs/figures/F11_station_timeseries.png", p11, width = 9, height = 4, dpi = 200)

write.csv(v, "outputs/tables/T2_station_validation.csv", row.names = FALSE)

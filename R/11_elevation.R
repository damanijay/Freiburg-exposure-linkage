# ============================================================
# 11_elevation.R - Copernicus DEM, elevation as a confounder
# ============================================================
library(rstac); library(terra); library(sf); library(dplyr)
library(exactextractr); library(ggplot2)

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
hex  <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
kids <- st_read("data_cache/cohort.gpkg", quiet = TRUE)
bb   <- as.numeric(st_bbox(st_transform(st_buffer(st_union(hex), 500), 4326)))

# 1. Download Copernicus DEM 30 m (two tiles, joined virtually)
dem_items <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  stac_search(collections = "cop-dem-glo-30", bbox = bb) |>
  post_request() |> items_sign(sign_fn = sign_planetary_computer())
urls <- paste0("/vsicurl/", assets_url(dem_items, asset_names = "data"))
cat("DEM tiles:", length(urls), "\n")
dem <- if (length(urls) > 1) vrt(urls) else rast(urls)
dem <- crop(dem, ext(bb[1], bb[3], bb[2], bb[4]))
dem <- project(dem, "EPSG:25832", method = "bilinear", res = 30)   # height varies smoothly
writeRaster(dem, "data_cache/dem_30m.tif", overwrite = TRUE)
print(global(mask(dem, vect(city)), c("min", "mean", "max"), na.rm = TRUE))   # ~200 to ~1280 m

# 2. Elevation per hexagon and at each child's home
hex$elev  <- exact_extract(dem, hex, "mean", progress = FALSE)
kids$elev <- terra::extract(dem, vect(kids))[, 2]
cat("Children's homes, elevation range (m):", round(range(kids$elev)), "\n")

# 3. Summer surface temperature per hexagon (Jun-Aug, gap-filled)
summer <- read.csv("data_cache/lst_hex_filled.csv") |> mutate(date = as.Date(date)) |>
  filter(format(date, "%m") %in% c("06", "07", "08")) |>
  group_by(hex_id) |> summarise(lst_summer = mean(lst_filled), .groups = "drop")

d <- hex |> st_drop_geometry() |> left_join(summer, by = "hex_id") |>
  mutate(green10 = green * 10,          # per 10 percentage points of green
         elev100 = elev / 100)          # per 100 m of height
cat("Correlation green vs elevation:", round(cor(d$green, d$elev), 2), "\n")

# 4. The confounding test
m1 <- lm(lst_summer ~ green10, data = d)
m2 <- lm(lst_summer ~ green10 + elev100, data = d)
res <- data.frame(
  model = c("1: green only", "2: green + elevation"),
  green_per_10pp = round(c(coef(m1)["green10"], coef(m2)["green10"]), 2),
  elev_per_100m  = round(c(NA, coef(m2)["elev100"]), 2),
  R2 = round(c(summary(m1)$r.squared, summary(m2)$r.squared), 2))
print(res)
cat("Green coefficient shrinks by",
    round(100 * (1 - coef(m2)["green10"] / coef(m1)["green10"])), "% after adjusting for elevation\n")
write.csv(res, "outputs/tables/T3_elevation_confounding.csv", row.names = FALSE)

# 5. Map
p12 <- ggplot(hex) + geom_sf(aes(fill = elev), colour = NA) +
  geom_sf(data = kids, size = 0.3, colour = "white", alpha = 0.6) +
  geom_sf(data = city, fill = NA, colour = "black") +
  scale_fill_viridis_c(name = "Elevation (m)", option = "cividis") +
  labs(title = "Elevation per hexagon, Freiburg", subtitle = "Copernicus DEM GLO-30; white dots = children") +
  theme_void()
print(p12); ggsave("outputs/figures/F12_elevation_hex.png", p12, width = 7, height = 7, dpi = 200)

# 6. Save (elevation now available to every later step)
st_write(hex,  "data_cache/hex_green_blue.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(kids, "data_cache/cohort.gpkg",         delete_dsn = TRUE, quiet = TRUE)

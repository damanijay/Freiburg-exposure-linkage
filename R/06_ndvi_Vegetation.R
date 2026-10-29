library(rstac); library(terra); library(sf); library(dplyr)
library(exactextractr); library(ggplot2)

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
hex  <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
area <- st_buffer(st_union(hex), 500)                         # hexagons + 500 m margin, in metres
bb   <- as.numeric(st_bbox(st_transform(area, 4326)))

# ---------- Section 1: search summer 2023, clouds below 20% ----------
s2 <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  stac_search(collections = "sentinel-2-l2a", bbox = bb,
              datetime = "2023-06-01/2023-08-31", limit = 200) |>
  ext_query("eo:cloud_cover" < 20) |>
  post_request() |> items_fetch() |>
  items_sign(sign_fn = sign_planetary_computer())

scenes <- data.frame(
  i     = seq_along(s2$features),
  date  = as.Date(substr(sapply(s2$features, \(f) f$properties$datetime), 1, 10)),
  tile  = sapply(s2$features, \(f) f$properties$`s2:mgrs_tile`),
  cloud = round(sapply(s2$features, \(f) f$properties$`eo:cloud_cover`), 1)) |>
  arrange(cloud)
print(scenes)                       # Freiburg may sit on 2 tiles; that's fine
pick <- head(scenes$i, 10)          # the 10 clearest scenes

# ---------- Section 2: one function = one scene → cloud-masked NDVI ----------
tmpl <- rast(ext(vect(area)), resolution = 10, crs = "EPSG:25832")   # common 10 m grid for all scenes

ndvi_scene <- function(i) {
  f  <- s2$features[[i]]
  rd <- function(a) { r <- rast(paste0("/vsicurl/", f$assets[[a]]$href))
  crop(r, project(vect(area), crs(r))) }        # read only our area
  red <- rd("B04"); nir <- rd("B08")
  scl <- resample(rd("SCL"), red, method = "near")                 # 20 m codes → 10 m, never blend
  red[red == 0] <- NA; nir[nir == 0] <- NA                          # 0 = no data
  off <- if (scenes$date[scenes$i == i] >= as.Date("2022-01-25")) 1000 else 0
  red <- (red - off) / 10000; nir <- (nir - off) / 10000            # stored value → reflectance
  ndvi <- (nir - red) / (nir + red)
  ndvi[scl %in% c(0, 1, 3, 8, 9, 10, 11)] <- NA                     # drop cloud, shadow, snow, defects
  ndvi <- clamp(ndvi, -1, 1, values = FALSE)
  project(ndvi, tmpl, method = "bilinear")                          # onto the common grid
}

# ---------- Section 3: process the scenes (the slow part) ----------
layers <- list()
for (i in pick) {
  r <- try(ndvi_scene(i), silent = TRUE)
  lab <- paste(scenes$date[scenes$i == i], scenes$tile[scenes$i == i])
  if (inherits(r, "try-error")) { message(lab, "  FAILED") ; next }
  names(r) <- lab; layers[[lab]] <- r; message(lab, "  OK")
}
ndvi_stack <- rast(layers)

# ---------- Section 4: median composite + checks ----------
ndvi_med <- median(ndvi_stack, na.rm = TRUE)
n_obs    <- sum(!is.na(ndvi_stack))                                 # clear views per pixel
ndvi_med <- mask(ndvi_med, vect(area))
writeRaster(ndvi_med, "data_cache/ndvi_summer2023.tif", overwrite = TRUE)

print(global(ndvi_med, c("min", "mean", "max"), na.rm = TRUE))      # expect mean ~0.5–0.7
print(global(n_obs, c("min", "mean"), na.rm = TRUE))                # every pixel should have ≥ 1
plot(ndvi_med, main = "Sentinel-2 NDVI, summer 2023 median (10 m)"); plot(st_geometry(city), add = TRUE)

# ---------- Section 5: NDVI per hexagon + agreement with WorldCover ----------
hex$ndvi <- exact_extract(ndvi_med, hex, "mean", progress = FALSE)
cat("Correlation NDVI vs WorldCover green share:", round(cor(hex$ndvi, hex$green, use = "complete.obs"), 2), "\n")

ggplot(hex) + geom_sf(aes(fill = ndvi), colour = NA) +
  geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.6) +
  scale_fill_viridis_c(name = "Mean NDVI") +
  labs(title = "Vegetation per hexagon, Freiburg", subtitle = "Sentinel-2 L2A, Jun–Aug 2023 median") +
  theme_void()
ggsave("outputs/figures/F3_ndvi_hex.png", width = 7, height = 7, dpi = 200)

st_write(hex, "data_cache/hex_green_blue.gpkg", delete_dsn = TRUE, quiet = TRUE)   # save progress

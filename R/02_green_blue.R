library(rstac); library(terra); library(sf); library(dplyr)
library(exactextractr); library(osmdata); library(ggplot2)

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
hex  <- st_read("data_cache/hex_1km.gpkg", quiet = TRUE)

# ---------- Part A: how much of each hexagon is inside Freiburg ----------
hex$hex_area <- as.numeric(st_area(hex))
inter <- st_intersection(hex["hex_id"], st_geometry(city))     # the inside part of each hexagon
inter$inside_area <- as.numeric(st_area(inter))
hex <- hex |>
  left_join(st_drop_geometry(inter)[, c("hex_id", "inside_area")], by = "hex_id") |>
  mutate(inside_frac = inside_area / hex_area)

summary(hex$inside_frac)                                        # 1 = fully inside
cat("Sum of inside areas (km²):", round(sum(hex$inside_area) / 1e6, 1), "\n")  # should be 152.9
cat("Total hexagon area (km²):", round(sum(hex$hex_area) / 1e6, 1), "\n")      # more than 152.9



# ---------- Part B: download ESA WorldCover 2021 for the hexagon area ----------
# Use the hexagons' extent (not just the city), so edge hexagons get full land cover
bb <- as.numeric(st_bbox(st_transform(st_buffer(st_union(hex), 500), 4326)))

wc_items <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  stac_search(collections = "esa-worldcover", bbox = bb) |>
  post_request() |> items_sign(sign_fn = sign_planetary_computer())
ids <- sapply(wc_items$features, \(f) f$id); print(ids)       # contains 2020 and 2021 versions
wc_items$features <- wc_items$features[grepl("2021", ids)]     # keep 2021

urls <- paste0("/vsicurl/", assets_url(wc_items, asset_names = "map"))
wc <- if (length(urls) > 1) vrt(urls) else rast(urls)           # join tiles if there are several
wc <- crop(wc, ext(bb[1], bb[3], bb[2], bb[4]))                 # read only our area
wc <- project(wc, "EPSG:25832", method = "near", res = 10)      # "near": never blend classes
writeRaster(wc, "data_cache/worldcover_2021.tif", overwrite = TRUE)

freq(wc)                                                        # how many 10 m cells per class
plot(wc, main = "ESA WorldCover 2021 (class codes)"); plot(st_geometry(city), add = TRUE)



# ---------- Part C: share of each class in every hexagon ----------
fr <- exact_extract(wc, hex, fun = "frac", append_cols = "hex_id", progress = FALSE)
share <- \(df, k) rowSums(df[, intersect(paste0("frac_", k), names(df)), drop = FALSE])

hex <- hex |> left_join(
  data.frame(hex_id   = fr$hex_id,
             green    = share(fr, c(10, 20, 30, 90)),   # trees, shrubs, grass, wetland
             cropland = share(fr, 40),
             built    = share(fr, 50),
             blue     = share(fr, 80)),
  by = "hex_id")

# ---------- Part C2: cap shares at 0-1 (removes floating-point rounding, e.g. 1.0000000000002) ----------
hex <- hex |> mutate(across(c(green, cropland, built, blue), \(x) pmin(pmax(x, 0), 1)))
stopifnot(all(hex$green >= 0 & hex$green <= 1))   # pipeline check: stops loudly if shares ever leave 0-1




# ---------- Part D: city-wide green share, three ways ----------
direct_fr <- exact_extract(wc, city, fun = "frac", progress = FALSE)
compare <- data.frame(
  method = c("Unweighted (every hexagon equal)", "Weighted by inside_frac", "Direct from city outline (truth)"),
  green  = round(c(mean(hex$green),
                   weighted.mean(hex$green, hex$inside_frac),
                   share(direct_fr, c(10, 20, 30, 90))), 3),
  built  = round(c(mean(hex$built),
                   weighted.mean(hex$built, hex$inside_frac),
                   share(direct_fr, 50)), 3))
print(compare)




# ---------- Part E: map ----------
ggplot(hex) +
  geom_sf(aes(fill = green), colour = NA) +
  geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.6) +
  scale_fill_viridis_c(name = "Green share", limits = c(0, 1)) +
  labs(title = "Green space per 1 km² hexagon, Freiburg",
       subtitle = "ESA WorldCover 2021: tree, shrub, grass and wetland (cropland excluded)") +
  theme_void()
ggsave("outputs/figures/F2_green_share_hex.png", width = 7, height = 7, dpi = 200)

st_write(hex, "data_cache/hex_green_blue.gpkg", delete_dsn = TRUE, quiet = TRUE)   # save progress





# ---------- Part F: OpenStreetMap as a second, independent source ----------
library(osmdata)
set_overpass_url("https://overpass-api.de/api/interpreter")   # main Overpass server
dir.create("data_cache/osm", showWarnings = FALSE)

# Download one feature type, retrying up to 3 times, and cache the result on disk
get_osm <- function(feature, tries = 3) {
  file <- file.path("data_cache/osm", paste0(gsub('[^a-z]', '_', feature), ".rds"))
  if (file.exists(file)) return(readRDS(file))                     # already downloaded: reuse it
  for (i in 1:tries) {
    res <- try(opq(bb, timeout = 300) |> add_osm_features(features = feature) |> osmdata_sf(),
               silent = TRUE)
    if (!inherits(res, "try-error")) { saveRDS(res, file); message(feature, ": OK"); return(res) }
    message(feature, ": try ", i, " failed, waiting 30 s")
    Sys.sleep(30)
  }
  stop("All tries failed for ", feature, ". Run the line again in a few minutes.")
}

green_feats <- c('"leisure"="park"', '"landuse"="forest"', '"natural"="wood"',
                 '"landuse"="meadow"', '"landuse"="grass"')
blue_feats  <- c('"natural"="water"', '"waterway"="riverbank"')

green_list <- lapply(green_feats, get_osm)
blue_list  <- lapply(blue_feats,  get_osm)

# Combine all polygons from several downloads, clean them, dissolve overlaps
to_poly <- function(lst) {
  parts <- unlist(lapply(lst, \(x) list(x$osm_polygons, x$osm_multipolygons)), recursive = FALSE)
  parts <- Filter(\(p) !is.null(p) && nrow(p) > 0, parts)
  g <- do.call(c, lapply(parts, st_geometry))
  st_sf(geometry = g) |> st_make_valid() |> st_transform(25832) |> st_union() |> st_sf()
}
green_osm <- to_poly(green_list)
blue_osm  <- to_poly(blue_list)

osm_share <- function(poly) {
  i <- st_intersection(hex["hex_id"], st_geometry(poly))
  out <- data.frame(hex_id = i$hex_id, a = as.numeric(st_area(i))) |>
    group_by(hex_id) |> summarise(a = sum(a), .groups = "drop")
  hex["hex_id"] |> st_drop_geometry() |> left_join(out, by = "hex_id") |>
    mutate(a = coalesce(a, 0)) |> pull(a)
}
hex$green_osm <- pmin(osm_share(green_osm) / hex$hex_area, 1)
hex$blue_osm  <- pmin(osm_share(blue_osm)  / hex$hex_area, 1)

cat("Agreement WorldCover vs OSM green (correlation):",
    round(cor(hex$green, hex$green_osm), 2), "\n")

st_write(hex, "data_cache/hex_green_blue.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(green_osm, "data_cache/osm_green.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(blue_osm,  "data_cache/osm_blue.gpkg",  delete_dsn = TRUE, quiet = TRUE)

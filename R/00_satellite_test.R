library(rstac); library(terra); library(sf)

# Load the city outline from Block 1
city <- st_read("data_cache/city.gpkg", quiet = TRUE)

# STAC searches in longitude/latitude, so convert the outline's bounding box to degrees (EPSG:4326)
bbox <- as.numeric(st_bbox(st_transform(city, 4326)))
bbox                                   # 4 numbers: west, south, east, north (around 7.6, 47.9, 8.0, 48.1)

# Search the catalogue: which collection, where, when
items <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  stac_search(collections = "modis-11A2-061", bbox = bbox,
              datetime = "2023-07-01/2023-07-31", limit = 50) |>
  post_request() |>
  items_sign(sign_fn = sign_planetary_computer())   # add the temporary access token

# The collection holds both Terra (MOD...) and Aqua (MYD...) images. Keep Terra only.
ids <- sapply(items$features, \(f) f$id)
print(head(ids))
items$features <- items$features[startsWith(ids, "MOD")]
cat("Terra composites in July 2023:", length(items$features), "\n")   # expect about 4

# Take the first composite and open its daytime temperature layer from the cloud
url <- assets_url(items, asset_names = "LST_Day_1km")[1]
r   <- rast(paste0("/vsicurl/", url))
r                                      # prints the info: huge tile, sinusoidal projection

# Cut out only Freiburg (the outline must be in the satellite's projection to match)
lst <- crop(r, vect(st_transform(city, crs(r))))

# Convert the stored numbers into °C
lst[lst == 0] <- NA                                            # 0 = no data
if (global(lst, "max", na.rm = TRUE)[1, 1] > 1000) lst <- lst * 0.02   # stored values → Kelvin
lst_c <- lst - 273.15                                          # Kelvin → °C

# Put it into the same metre-based system as everything else
lst_c <- project(lst_c, "EPSG:25832")

# Check the values, then draw the map
print(global(lst_c, c("min", "mean", "max"), na.rm = TRUE))   # expect roughly 25–45 °C
plot(lst_c, main = "Freiburg daytime surface temperature, July 2023 (°C)")
plot(st_geometry(city), add = TRUE)

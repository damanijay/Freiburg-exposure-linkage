install.packages("renv")
renv::consent(provided = TRUE)
renv::init()

pkgs <- c("sf","terra","rstac","exactextractr","osmdata","dplyr","tidyr",
          "purrr","zoo","ggplot2","lubridate","here","targets","tarchetypes",
          "DBI","duckdb","testthat","sandwich","lmtest","rmarkdown")
renv::install(pkgs)
renv::snapshot()

# Create folder structure
for (d in c("R","data_cache","outputs/figures","outputs/tables","db","docs","tests"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Check spatial libraries
sf::sf_extSoftVersion()
terra::gdal()


library(sf); library(dplyr)
crs_m <- 25832                      # metres-based CRS for SW Germany

# 1. Ask OpenStreetMap's geocoder for Freiburg's boundary polygon
city <- osmdata::getbb("Freiburg im Breisgau, Germany", format_out = "sf_polygon")
if (!inherits(city, "sf")) city <- do.call(rbind, city)

# 2. Clean it: fix invalid geometry, merge pieces into one, convert to metres
city <- city |> st_make_valid() |> st_union() |> st_as_sf() |> st_transform(crs_m)

# 3. Hexagon grid. A hexagon 1,075 m across (flat side to flat side) has an area of about 1 km²,
#    the same as one MODIS pixel. Hexagons have equal distance to all neighbours,
#    which is why spatial epidemiology prefers them to squares.
hex <- st_make_grid(city, cellsize = 1075, square = FALSE) |>
  st_as_sf() |> st_filter(city) |> mutate(hex_id = row_number())

# 4. Quality checks: always check before you move on
cat("City area (km²):", round(as.numeric(st_area(city)) / 1e6, 1), "\n")  # expect ~153
cat("Hexagons:", nrow(hex), "\n")                                         # expect ~150-190
stopifnot(st_crs(city)$epsg == crs_m)

st_write(city, "data_cache/city.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(hex,  "data_cache/hex_1km.gpkg", delete_dsn = TRUE, quiet = TRUE)
plot(st_geometry(hex), border = "grey"); plot(st_geometry(city), add = TRUE, lwd = 2)

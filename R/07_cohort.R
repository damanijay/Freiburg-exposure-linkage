# ============================================================
# 07_cohort.R  -  Simulated cohort of 500 children, Freiburg
# ============================================================
library(terra); library(sf); library(dplyr); library(ggplot2)
set.seed(2026)                                   # same children every run

city <- st_read("data_cache/city.gpkg", quiet = TRUE)
hex  <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
wc   <- rast("data_cache/worldcover_2021.tif")

# 1. Sampling frame: built-up cells (class 50) inside the city, everything else NA
bu <- classify(wc, cbind(50, 1), others = NA) |> mask(vect(city))

# 2. Draw homes: sample generously, check there are enough, keep 500
pts <- spatSample(bu, size = 5000, method = "random", na.rm = TRUE, as.points = TRUE) |>
  st_as_sf()
cat("Points landing on built-up cells:", nrow(pts), "\n")
stopifnot(nrow(pts) >= 500)                      # fail loudly if the frame is too small

kids <- pts |> slice_head(n = 500) |> select(geometry) |>
  mutate(pt_id = row_number())                   # temporary ID per home

# 3. Spatial join: which hexagon each child lives in
#    (a home on a shared edge can match 2 hexagons, so keep one row per home)
kids <- st_join(kids, hex[c("hex_id", "built")]) |>
  distinct(pt_id, .keep_all = TRUE) |>
  select(-pt_id)

# 4. Child characteristics
kids <- kids |> mutate(
  child_id  = sprintf("C%03d", row_number()),
  exam_date = as.Date("2022-01-01") + sample(0:729, n(), replace = TRUE),  # 2022-2023
  age       = round(runif(n(), 2, 9.9), 1),                                # I.Family age range
  sex       = sample(c("f", "m"), n(), replace = TRUE),
  # Planted truth: SES falls 0.5 SD per 1 SD more built-up, plus individual noise
  ses_z     = as.numeric(scale(-0.5 * as.numeric(scale(built)) + rnorm(n(), 0, 1)))
) |> select(child_id, hex_id, exam_date, age, sex, ses_z, built)

# 5. Checks
stopifnot(nrow(kids) == 500, !anyNA(kids$hex_id), !anyDuplicated(kids$child_id))
print(summary(kids$exam_date))                   # 2022-01-01 to 2023-12-31
cat("Children per quarter:\n"); print(table(quarters(kids$exam_date)))
cat("Hexagons with at least one child:", n_distinct(kids$hex_id), "\n")
cat("Check planted SES link (should be about -0.4 to -0.5):",
    round(cor(kids$ses_z, kids$built), 2), "\n")

# 6. Map
p <- ggplot() +
  geom_sf(data = hex, fill = "grey95", colour = "grey80") +
  geom_sf(data = kids, aes(colour = ses_z), size = 0.8) +
  geom_sf(data = city, fill = NA, colour = "black") +
  scale_colour_gradient2(name = "SES (z)", low = "firebrick", mid = "grey70", high = "steelblue") +
  labs(title = "Simulated cohort: 500 children, Freiburg",
       subtitle = "Homes sampled on built-up cells (ESA WorldCover 2021)") +
  theme_void()
print(p)
ggsave("outputs/figures/F4_cohort.png", p, width = 7, height = 7, dpi = 200)

# 7. Save
st_write(kids, "data_cache/cohort.gpkg", delete_dsn = TRUE, quiet = TRUE)
message("Cohort saved: ", nrow(kids), " children")

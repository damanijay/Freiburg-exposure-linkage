# ============================================================
# 08_lst_cube.R  -  LST data cube, gap diagnosis, child windows
# ============================================================
library(terra); library(sf); library(dplyr); library(tidyr)
library(exactextractr); library(ggplot2)

hex  <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
kids <- st_read("data_cache/cohort.gpkg", quiet = TRUE)
inv  <- read.csv("data_cache/lst_inventory.csv") |> mutate(date = as.Date(date))

# ---------- Section 1: build the hexagon x week cube ----------
files <- list.files("data_cache/lst", pattern = "\\.tif$", full.names = TRUE)
dates <- as.Date(sub(".*lst_(.*)\\.tif", "\\1", files))
lst   <- rast(files)                              # 130 layers, one per week
valid <- 1 - is.na(lst)                           # 1 = valid pixel, 0 = cloud / no data

hex_s <- st_transform(hex["hex_id"], crs(lst))    # hexagons into the satellite's projection
m_lst <- exact_extract(lst,   hex_s, "mean", progress = FALSE)   # mean LST of valid pixels
m_val <- exact_extract(valid, hex_s, "mean", progress = FALSE)   # share of hexagon that is valid
colnames(m_lst) <- colnames(m_val) <- as.character(dates)

to_long <- \(m, v) data.frame(hex_id = hex$hex_id, m, check.names = FALSE) |>
  pivot_longer(-hex_id, names_to = "date", values_to = v) |> mutate(date = as.Date(date))

cube <- expand_grid(hex_id = hex$hex_id, date = inv$date) |>       # all 138 weeks, incl. missing
  left_join(to_long(m_lst, "lst"),        by = c("hex_id", "date")) |>
  left_join(to_long(m_val, "valid_frac"), by = c("hex_id", "date")) |>
  mutate(valid_frac = coalesce(valid_frac, 0),
         lst    = ifelse(valid_frac >= 0.5, lst, NA),               # the 50% rule
         slot   = (as.integer(format(date, "%j")) - 1) %/% 8 + 1,  # week-of-year slot 1-46
         season = ifelse(as.integer(format(date, "%m")) %in% c(10:12, 1:3),
                         "cold (Oct-Mar)", "warm (Apr-Sep)"))

stopifnot(nrow(cube) == nrow(hex) * nrow(inv))
cat("Hexagon-weeks:", nrow(cube), " observed:", sum(!is.na(cube$lst)),
    " missing:", round(100 * mean(is.na(cube$lst)), 1), "%\n")
write.csv(cube, "data_cache/lst_hex_cube.csv", row.names = FALSE)

# ---------- Section 2: when and where is data missing? ----------
gap_month <- cube |> mutate(month = as.integer(format(date, "%m"))) |>
  group_by(month) |> summarise(pct_missing = round(100 * mean(is.na(lst)), 1))
print(gap_month)

p5 <- ggplot(gap_month, aes(factor(month), pct_missing)) +
  geom_col(fill = "steelblue") +
  labs(x = "Month", y = "% hexagon-weeks missing",
       title = "MODIS LST gaps by month, Freiburg 2021-2023") + theme_minimal()
print(p5); ggsave("outputs/figures/F5_gaps_by_month.png", p5, width = 7, height = 4, dpi = 200)

gap_hex <- cube |> filter(season == "cold (Oct-Mar)") |>
  group_by(hex_id) |> summarise(pct_missing_cold = 100 * mean(is.na(lst)))
hex_gap <- hex |> left_join(gap_hex, by = "hex_id") |>
  mutate(has_children = hex_id %in% kids$hex_id)

p6 <- ggplot(hex_gap) + geom_sf(aes(fill = pct_missing_cold), colour = NA) +
  scale_fill_viridis_c(name = "% missing\n(Oct-Mar)", option = "magma") +
  labs(title = "Where winter satellite gaps fall") + theme_void()
print(p6); ggsave("outputs/figures/F6_gaps_cold_map.png", p6, width = 7, height = 7, dpi = 200)

hex_gap |> st_drop_geometry() |> group_by(has_children) |>
  summarise(mean_pct_missing_cold = round(mean(pct_missing_cold), 1)) |> print()

# ---------- Section 3: each child's 7-, 30-, 365-day window ----------
periods <- inv |> transmute(date, start = date,
                            end = pmin(date + 7, as.Date(paste0(format(date, "%Y"), "-12-31"))))   # last block stops on 31 Dec

kid_tab <- kids |> st_drop_geometry() |>
  transmute(child_id, hex_id, exam_date = as.Date(exam_date))

link <- crossing(kid_tab, window = c(7, 30, 365)) |>
  mutate(wstart = exam_date - window, wend = exam_date - 1) |>
  cross_join(periods) |>
  mutate(overlap = as.numeric(pmin(end, wend) - pmax(start, wstart)) + 1) |>
  filter(overlap > 0) |>
  left_join(cube |> select(hex_id, date, lst), by = c("hex_id", "date"))

chk <- link |> group_by(child_id, window) |> summarise(days = sum(overlap), .groups = "drop")
stopifnot(all(chk$days == chk$window))           # overlap days must add up to the window

expo_naive <- link |> group_by(child_id, window, exam_date) |>
  summarise(coverage  = sum(overlap[!is.na(lst)]) / first(window),
            lst_naive = if (any(!is.na(lst))) weighted.mean(lst[!is.na(lst)], overlap[!is.na(lst)])
            else NA_real_,
            .groups = "drop") |>
  mutate(quarter = quarters(exam_date))

expo_naive |> group_by(window) |>
  summarise(pct_no_data      = round(100 * mean(is.na(lst_naive)), 1),
            pct_coverage_lt50 = round(100 * mean(coverage < 0.5), 1),
            mean_coverage     = round(mean(coverage), 2)) |> print()

expo_naive |> group_by(window, quarter) |>
  summarise(pct_coverage_lt50 = round(100 * mean(coverage < 0.5), 1), .groups = "drop") |>
  pivot_wider(names_from = quarter, values_from = pct_coverage_lt50) |> print()

write.csv(link,       "data_cache/window_link.csv",    row.names = FALSE)
write.csv(expo_naive, "data_cache/exposure_naive.csv", row.names = FALSE)

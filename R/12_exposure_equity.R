# ============================================================
# 12_exposure_equity.R - child exposures, when vs where, equity
# ============================================================
library(sf); library(terra); library(dplyr); library(tidyr)
library(exactextractr); library(ggplot2); library(sandwich); library(lmtest)

kids      <- st_read("data_cache/cohort.gpkg", quiet = TRUE)
wc        <- rast("data_cache/worldcover_2021.tif")
ndvi      <- rast("data_cache/ndvi_summer2023.tif")
green_osm <- st_read("data_cache/osm_green.gpkg", quiet = TRUE)
expo      <- read.csv("data_cache/exposure_final.csv")

# ---------- Section 1: green exposure in a 300 m home buffer ----------
buf <- st_buffer(kids, 300)
fr  <- exact_extract(wc, buf, "frac", progress = FALSE)
share <- \(df, k) rowSums(df[, intersect(paste0("frac_", k), names(df)), drop = FALSE])
kids$green300 <- pmin(share(fr, c(10, 20, 30, 90)), 1)      # trees, shrubs, grass, wetland
kids$ndvi300  <- exact_extract(ndvi, buf, "mean", progress = FALSE)

# Park access: public green >= 0.5 ha within 300 m (WHO indicator)
g <- st_geometry(green_osm)
g <- suppressWarnings(st_collection_extract(g, "POLYGON"))
parks <- st_sf(geometry = st_cast(st_cast(g, "MULTIPOLYGON"), "POLYGON"))   # split back into areas
parks <- parks[as.numeric(st_area(parks)) >= 5000, ]                         # >= 0.5 ha
nn <- st_nearest_feature(kids, parks)
kids$dist_park   <- as.numeric(st_distance(kids, parks[nn, ], by_element = TRUE))
kids$park_access <- as.numeric(kids$dist_park <= 300)
cat("Green areas >= 0.5 ha:", nrow(parks), "\n")
cat("Children with park access (<= 300 m):", round(100 * mean(kids$park_access), 1), "%\n")
summary(kids[c("green300", "ndvi300", "dist_park")] |> st_drop_geometry())

# ---------- Section 2: one row per child (the final exposure table) ----------
lst_wide <- expo |> select(child_id, window, lst_filled, coverage) |>
  pivot_wider(names_from = window, values_from = c(lst_filled, coverage),
              names_glue = "{.value}_{window}")            # e.g. lst_filled_365, coverage_7

final <- kids |> st_drop_geometry() |> left_join(lst_wide, by = "child_id")
stopifnot(nrow(final) == 500, !anyNA(final$lst_filled_365), !anyNA(final$green300))
write.csv(final, "data_cache/exposure_child_final.csv", row.names = FALSE)
cat("Final exposure table:", nrow(final), "children x", ncol(final), "columns\n")

# ---------- Section 3: when vs where ----------
ww <- expo |> left_join(final |> select(child_id, hex_id), by = "child_id") |>
  mutate(exam_month = format(as.Date(exam_date), "%Y-%m"))

when_where <- ww |> group_by(window) |>
  group_modify(\(d, k) data.frame(
    R2_time  = summary(lm(lst_filled ~ factor(exam_month), data = d))$adj.r.squared,  # WHEN
    R2_place = summary(lm(lst_filled ~ factor(hex_id),    data = d))$adj.r.squared)) |>  # WHERE
  mutate(across(starts_with("R2"), \(x) round(x, 2)))
print(when_where)
write.csv(when_where, "outputs/tables/T4_when_vs_where.csv", row.names = FALSE)

cat("Correlation between windows:\n")
print(round(cor(final[c("lst_filled_7", "lst_filled_30", "lst_filled_365")]), 2))

p13 <- ww |> mutate(window = factor(paste(window, "days"), levels = c("7 days", "30 days", "365 days"))) |>
  ggplot(aes(as.Date(exam_date), lst_filled)) + geom_point(alpha = 0.4, size = 0.8) +
  facet_wrap(~window, ncol = 1) +
  labs(title = "Same children, three exposure windows",
       subtitle = "Short windows track the exam season; the 365-day window tracks the neighbourhood",
       x = "Exam date", y = "LST exposure (°C)") + theme_minimal()
print(p13); ggsave("outputs/figures/F13_windows_when_where.png", p13, width = 8, height = 7, dpi = 200)

# ---------- Section 4a: does the method recover the planted SES link? ----------
final <- final |> mutate(built_z = as.numeric(scale(built)))
m0 <- lm(ses_z ~ built_z, data = final)
V_cl <- vcovCL(m0, cluster = ~hex_id)                       # children clustered in 80 hexagons
ci_naive <- confint(m0)["built_z", ]
ci_clust <- coefci(m0, vcov. = V_cl)["built_z", ]
truth <- -0.5 / sqrt(1.25)                                  # planted slope after rescaling

recovery <- data.frame(
  estimate = round(coef(m0)["built_z"], 3),
  truth    = round(truth, 3),
  SE_naive = round(sqrt(vcov(m0)["built_z", "built_z"]), 3),
  SE_clust = round(sqrt(V_cl["built_z", "built_z"]), 3),
  CI_naive = paste(round(ci_naive, 3), collapse = " to "),
  CI_clust = paste(round(ci_clust, 3), collapse = " to "),
  truth_in_naive_CI = truth >= ci_naive[1] & truth <= ci_naive[2],
  truth_in_clust_CI = truth >= ci_clust[1] & truth <= ci_clust[2])
print(recovery)
write.csv(recovery, "outputs/tables/T5_truth_recovery.csv", row.names = FALSE)

# ---------- Section 4b: equity - do lower-SES children get more heat, less green? ----------
outs <- c(lst_filled_365 = "Annual LST (°C)",
          green300       = "Green share, 300 m",
          ndvi300        = "NDVI, 300 m",
          park_access    = "Park >= 0.5 ha within 300 m (prob.)")

equity <- bind_rows(lapply(names(outs), \(y) {
  m  <- lm(reformulate("ses_z", y), data = final)          # unadjusted on purpose (descriptive)
  V  <- vcovCL(m, cluster = ~hex_id)
  ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V)
  data.frame(outcome = outs[[y]],
             per_1SD_higher_SES = signif(ct["ses_z", 1], 3),
             CI_low  = signif(ci["ses_z", 1], 3),
             CI_high = signif(ci["ses_z", 2], 3),
             p_clustered = signif(ct["ses_z", 4], 2))
}))
print(equity)
write.csv(equity, "outputs/tables/T6_equity.csv", row.names = FALSE)

final$ses_group <- factor(ntile(final$ses_z, 3), labels = c("Low SES", "Middle", "High SES"))
final |> group_by(ses_group) |>
  summarise(annual_LST = round(mean(lst_filled_365), 2),
            green300   = round(mean(green300), 2),
            park_pct   = round(100 * mean(park_access), 1), n = n()) |> print()

p14 <- final |> select(ses_group, `Annual LST (°C)` = lst_filled_365, `Green share 300 m` = green300) |>
  pivot_longer(-ses_group) |>
  ggplot(aes(ses_group, value, fill = ses_group)) + geom_boxplot(alpha = 0.7) +
  facet_wrap(~name, scales = "free_y") +
  labs(title = "Environmental exposure by SES group (simulated cohort, real exposure data)",
       x = NULL, y = NULL) + theme_minimal() + theme(legend.position = "none")
print(p14); ggsave("outputs/figures/F14_equity.png", p14, width = 8, height = 4.5, dpi = 200)

# ============================================================
# 09_interpolation.R  -  v1.1
# Outlier screen + robust anomaly, validate gap-filling, fill, re-link
# ============================================================
library(dplyr); library(tidyr); library(zoo); library(ggplot2)
set.seed(2026)

cube <- read.csv("data_cache/lst_hex_cube.csv") |> mutate(date = as.Date(date))

# ---------- Section 1: matrix form + climatology function ----------
M <- cube |> select(hex_id, date, lst) |>
  pivot_wider(names_from = hex_id, values_from = lst) |> arrange(date)
X_raw <- as.matrix(M[, -1])                                 # rows = 138 weeks, cols = hexagons
slot  <- (as.integer(format(M$date, "%j")) - 1) %/% 8 + 1  # week-of-year slot 1-46

clim_of <- function(Z) {                                    # usual value: same slot +-1, all years
  C <- matrix(NA_real_, nrow(Z), ncol(Z))
  for (k in 1:46) {
    nb <- ((k - 2):k) %% 46 + 1                             # slots k-1, k, k+1 (wraps round the year)
    cm <- colMeans(Z[slot %in% nb, , drop = FALSE], na.rm = TRUE)
    C[slot == k, ] <- matrix(cm, sum(slot == k), ncol(Z), byrow = TRUE)
  }
  C[is.nan(C)] <- NA; C
}

# ---------- Section 2: spatial outlier screen (NEW in v1.1) ----------
C0 <- clim_of(X_raw)
A0 <- X_raw - C0                                    # how unusual each value is for ITS hexagon
D  <- A0 - apply(A0, 1, median, na.rm = TRUE)       # ...compared with the city that same week
sigma <- mad(D, na.rm = TRUE)                       # robust standard deviation of D
flag  <- !is.na(D) & abs(D) > 4 * sigma             # > 4 robust SD = almost surely an error

cat("Robust SD of deviations:", round(sigma, 2), "°C\n")
cat("Flagged as outliers:", sum(flag), "of", sum(!is.na(X_raw)), "observed values (",
    round(100 * sum(flag) / sum(!is.na(X_raw)), 2), "% )\n")

outliers <- data.frame(date   = M$date[row(X_raw)[flag]],
                       hex_id = as.integer(colnames(X_raw)[col(X_raw)[flag]]),
                       lst    = round(X_raw[flag], 1),
                       deviation = round(D[flag], 1))
outliers |> count(date, name = "n_hexagons") |> arrange(desc(n_hexagons)) |> head(10) |> print()
cat("Station hexagon 107 on 2023-11-01 flagged:",
    any(outliers$hex_id == 107 & outliers$date == as.Date("2023-11-01")), "\n")
write.csv(outliers, "outputs/tables/T1b_outliers_removed.csv", row.names = FALSE)

X <- X_raw; X[flag] <- NA                           # cleaned matrix: used for everything below

# ---------- Section 3: the three gap-filling methods ----------
fill_clim   <- function(Z) { C <- clim_of(Z); Z[is.na(Z)] <- C[is.na(Z)]; Z }
fill_linear <- function(Z) apply(Z, 2, \(x) if (sum(!is.na(x)) >= 2) na.approx(x, rule = 2) else x)
fill_anom   <- function(Z, min_frac = 0.2) {            # v1.1: robust anomaly
  C <- clim_of(Z)
  A <- Z - C
  a <- apply(A, 1, median, na.rm = TRUE)                # median: one bad hexagon can't drag it
  frac <- rowMeans(!is.na(A))                           # share of hexagons observed that week
  a[frac < min_frac | is.na(a)] <- NA                   # too few: treat the week as fully missing
  a <- na.approx(a, rule = 2)                           # borrow the anomaly from neighbouring weeks
  P <- C + a
  Z[is.na(Z)] <- P[is.na(Z)]; Z
}
methods <- list(climatology = fill_clim, linear = fill_linear, clim_anomaly = fill_anom)

# ---------- Section 4: hide known values, refill, score ----------
obs_idx <- which(!is.na(X))
mask_random <- function() sample(obs_idx, round(0.1 * length(obs_idx)))
mask_weeks  <- function() {
  full <- which(rowSums(!is.na(X)) > 0.9 * ncol(X))         # weeks that are almost complete
  w    <- sample(full, round(0.1 * nrow(X)))
  idx  <- as.vector(outer(w, (0:(ncol(X) - 1)) * nrow(X), "+"))   # every hexagon in those weeks
  idx[!is.na(X[idx])]
}

res <- list()
for (scheme in c("random_pixels", "whole_weeks")) for (r in 1:20) {
  idx <- if (scheme == "random_pixels") mask_random() else mask_weeks()
  Z <- X; Z[idx] <- NA                                      # methods only ever see Z (no leakage)
  for (m in names(methods)) {
    F <- methods[[m]](Z)
    res[[length(res) + 1]] <- data.frame(scheme, rep = r, method = m,
                                         row = (idx - 1) %% nrow(X) + 1, truth = X[idx], pred = F[idx])
  }
}
val <- bind_rows(res) |>
  mutate(season = ifelse(as.integer(format(M$date[row], "%m")) %in% c(10:12, 1:3),
                         "cold (Oct-Mar)", "warm (Apr-Sep)"))

scores <- val |> group_by(scheme, season, method) |>
  summarise(RMSE = round(sqrt(mean((pred - truth)^2, na.rm = TRUE)), 2),
            MAE  = round(mean(abs(pred - truth), na.rm = TRUE), 2),
            bias = round(mean(pred - truth, na.rm = TRUE), 2), .groups = "drop") |>
  arrange(scheme, season, RMSE)
print(scores, n = Inf)
write.csv(scores, "outputs/tables/T1_interpolation_validation_v1.1.csv", row.names = FALSE)

p7 <- ggplot(scores, aes(method, RMSE, fill = method)) + geom_col() +
  facet_grid(season ~ scheme) +
  labs(title = "Gap-filling validation v1.1: error on hidden real values (lower is better)",
       y = "RMSE (°C)", x = NULL) + theme_minimal() + theme(legend.position = "none")
print(p7); ggsave("outputs/figures/F7_interpolation_rmse.png", p7, width = 8, height = 5, dpi = 200)

# Pre-specified rule (unchanged): lowest cold-season RMSE averaged over both schemes
best <- scores |> filter(season == "cold (Oct-Mar)") |> group_by(method) |>
  summarise(RMSE = mean(RMSE)) |> slice_min(RMSE, n = 1) |> pull(method)
cat("Selected method:", best, "\n")

# ---------- Section 5: fill the real gaps, flag provenance ----------
F_all <- methods[[best]](X)
src <- ifelse(flag, "removed_outlier",
              ifelse(is.na(X_raw), paste0("filled_", best), "observed"))
colnames(src) <- colnames(X_raw)

to_long <- \(mat, v) data.frame(date = M$date, mat, check.names = FALSE) |>
  pivot_longer(-date, names_to = "hex_id", values_to = v) |> mutate(hex_id = as.integer(hex_id))

filled <- to_long(F_all, "lst_filled") |>
  left_join(to_long(src, "source"), by = c("date", "hex_id")) |>
  left_join(cube |> select(hex_id, date, lst), by = c("hex_id", "date"))   # lst = raw value
print(table(filled$source))
stopifnot(!anyNA(filled$lst_filled))                        # no gaps left
write.csv(filled, "data_cache/lst_hex_filled.csv", row.names = FALSE)

# One hexagon through 2022: the October hole being filled
h <- as.integer(names(sort(table(read.csv("data_cache/window_link.csv")$hex_id), decreasing = TRUE))[1])
p9 <- filled |> filter(hex_id == h, format(date, "%Y") == "2022") |>
  ggplot(aes(date, lst_filled, colour = source)) + geom_line(colour = "grey70") + geom_point(size = 2) +
  labs(title = paste("Hexagon", h, "(most children), 2022: observed vs filled (v1.1)"),
       y = "LST (°C)", x = NULL) + theme_minimal()
print(p9); ggsave("outputs/figures/F9_hex_timeseries.png", p9, width = 8, height = 4, dpi = 200)

# ---------- Section 6: re-link children with filled data ----------
expo_naive <- read.csv("data_cache/exposure_naive.csv")
expo_filled <- read.csv("data_cache/window_link.csv") |> mutate(date = as.Date(date)) |>
  left_join(filled |> select(hex_id, date, lst_filled), by = c("hex_id", "date")) |>
  group_by(child_id, window) |>
  summarise(lst_filled = weighted.mean(lst_filled, overlap), .groups = "drop")

cmp <- expo_naive |> left_join(expo_filled, by = c("child_id", "window")) |>
  mutate(diff = lst_filled - lst_naive)

cmp |> group_by(window) |>
  summarise(children_rescued = sum(is.na(lst_naive)),
            mean_diff        = round(mean(diff, na.rm = TRUE), 2),
            mean_abs_change  = round(mean(abs(diff), na.rm = TRUE), 2),
            max_abs_change   = round(max(abs(diff), na.rm = TRUE), 2)) |> print()

p8 <- cmp |> filter(window == 30) |>
  ggplot(aes(lst_naive, lst_filled, colour = coverage)) + geom_point(alpha = 0.7) +
  geom_abline(linetype = 2) + scale_colour_viridis_c() +
  labs(title = "30-day exposure: gaps ignored vs gaps filled (v1.1)",
       x = "Naive LST (°C, gaps ignored)", y = "Filled LST (°C)") + theme_minimal()
print(p8); ggsave("outputs/figures/F8_naive_vs_filled.png", p8, width = 6, height = 5, dpi = 200)

write.csv(cmp, "data_cache/exposure_final.csv", row.names = FALSE)

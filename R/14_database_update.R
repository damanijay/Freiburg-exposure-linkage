# ============================================================
# 14_database_update.R - incremental update: add new MODIS weeks (run 2)
# Same frozen v1.1 rules; earlier rows are never touched.
# ============================================================
library(DBI); library(duckdb); library(dplyr); library(tidyr)
library(terra); library(sf); library(exactextractr); library(zoo); library(testthat)

db_path <- "db/freiburg_exposure.duckdb"

# ---------- Section 1: tests BEFORE the update ----------
cat("\n===== TESTS BEFORE UPDATE =====\n")
test_file("tests/test_database.R")

con <- dbConnect(duckdb(), dbdir = db_path)
q <- function(sql) dbGetQuery(con, sql)

# ---------- Section 2: what is new in staging? ----------
files  <- list.files("data_cache/lst_incoming", pattern = "\\.tif$", full.names = TRUE)
fdates <- as.Date(sub(".*lst_(.*)\\.tif", "\\1", files))
existing <- q("SELECT date FROM lst_composite")$date
keep <- !fdates %in% existing
files <- files[keep]; fdates <- fdates[keep]
if (length(files) == 0) {
  dbDisconnect(con)
  stop("Nothing new in staging - the database is already up to date (this is expected on a rerun).")
}

# Expected 8-day calendar for the new year(s): restarts on 1 January every year
yrs <- sort(unique(format(fdates, "%Y")))
expected <- do.call(c, lapply(yrs, \(y) seq(as.Date(paste0(y, "-01-01")), by = 8, length.out = 46)))
new_dates <- expected[!expected %in% existing]
cat("New files in staging:", length(files), " | expected new weeks:", length(new_dates), "\n")

# ---------- Section 3: quality gate - refuse bad data before it enters ----------
ref <- rast(list.files("data_cache/lst", pattern = "\\.tif$", full.names = TRUE)[1])
same_grid <- sapply(files, \(f) isTRUE(compareGeom(ref, rast(f), stopOnError = FALSE)))
stopifnot("New files must use the same grid as the database" = all(same_grid),
          "New files must sit on the 8-day calendar"         = all(fdates %in% new_dates))
cat("Quality gate passed: same grid (", ncell(ref), "cells), dates on the calendar\n")

# ---------- Section 4: hexagon means for the new weeks (same 50% rule) ----------
hex   <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
r_new <- rast(files)
hex_s <- st_transform(hex["hex_id"], crs(r_new))
m_lst <- exact_extract(r_new,            hex_s, "mean", progress = FALSE)
m_val <- exact_extract(1 - is.na(r_new), hex_s, "mean", progress = FALSE)
colnames(m_lst) <- colnames(m_val) <- as.character(fdates)

to_long <- \(m, v) data.frame(hex_id = hex$hex_id, m, check.names = FALSE) |>
  pivot_longer(-hex_id, names_to = "date", values_to = v) |> mutate(date = as.Date(date))
new_long <- expand_grid(hex_id = hex$hex_id, date = new_dates) |>
  left_join(to_long(m_lst, "lst"),        by = c("hex_id", "date")) |>
  left_join(to_long(m_val, "valid_frac"), by = c("hex_id", "date")) |>
  mutate(valid_frac = coalesce(valid_frac, 0),
         lst = ifelse(valid_frac >= 0.5, lst, NA))

# ---------- Section 5: same outlier screen + gap-filling (frozen v1.1) ----------
old <- q("SELECT hex_id, date, lst_raw, source FROM lst_weekly")
old_n_before <- q("SELECT COUNT(*) AS n, SUM(lst_filled) AS s FROM lst_weekly WHERE run_id = 1")

all_v <- bind_rows(old |> transmute(hex_id, date, v = lst_raw),
                   new_long |> transmute(hex_id, date, v = lst)) |> arrange(date, hex_id)
all_f <- bind_rows(old |> transmute(hex_id, date, f = source == "removed_outlier"),
                   new_long |> transmute(hex_id, date, f = FALSE)) |> arrange(date, hex_id)
M  <- pivot_wider(all_v, names_from = hex_id, values_from = v)
Fl <- pivot_wider(all_f, names_from = hex_id, values_from = f)
stopifnot(identical(names(M), names(Fl)), identical(M$date, Fl$date))

X_raw    <- as.matrix(M[, -1]); old_flag <- as.matrix(Fl[, -1])
new_rows <- M$date %in% new_dates
slot     <- (as.integer(format(M$date, "%j")) - 1) %/% 8 + 1

clim_of <- function(Z) {
  C <- matrix(NA_real_, nrow(Z), ncol(Z))
  for (k in 1:46) {
    nb <- ((k - 2):k) %% 46 + 1
    cm <- colMeans(Z[slot %in% nb, , drop = FALSE], na.rm = TRUE)
    C[slot == k, ] <- matrix(cm, sum(slot == k), ncol(Z), byrow = TRUE)
  }
  C[is.nan(C)] <- NA; C
}
fill_anom <- function(Z, min_frac = 0.2) {
  C <- clim_of(Z); A <- Z - C
  a <- apply(A, 1, median, na.rm = TRUE)
  a[rowMeans(!is.na(A)) < min_frac | is.na(a)] <- NA
  a <- na.approx(a, rule = 2)
  P <- C + a; Z[is.na(Z)] <- P[is.na(Z)]; Z
}

# Outlier screen: same rule, threshold taken from the existing data (frozen)
C0 <- clim_of(X_raw); A0 <- X_raw - C0
D  <- A0 - apply(A0, 1, median, na.rm = TRUE)
sigma <- mad(D[!new_rows, ], na.rm = TRUE)
flag_new <- !is.na(D) & abs(D) > 4 * sigma
flag_new[!new_rows, ] <- FALSE
cat("Robust SD (frozen threshold):", round(sigma, 2), "°C | new outliers flagged:", sum(flag_new), "\n")

Xc <- X_raw; Xc[old_flag | flag_new] <- NA
Fm <- fill_anom(Xc)
src <- ifelse(flag_new, "removed_outlier", ifelse(is.na(X_raw), "filled_clim_anomaly", "observed"))

# ---------- Section 6: assemble the new rows ----------
run_id <- as.integer(q("SELECT MAX(run_id) AS m FROM run_log")$m) + 1L
nr <- which(new_rows); hex_ids <- as.integer(colnames(X_raw))
lst_new <- data.frame(
  hex_id = rep(hex_ids, each = length(nr)), date = rep(M$date[nr], times = ncol(X_raw)),
  lst_raw = as.vector(X_raw[nr, ]), lst_filled = as.vector(Fm[nr, ]),
  source = as.vector(src[nr, ]), method_version = "v1.1", run_id = run_id)
stopifnot(!anyNA(lst_new$lst_filled),
          all(is.na(lst_new$lst_raw) | (lst_new$lst_raw > -40 & lst_new$lst_raw < 70)))

fl <- which(flag_new, arr.ind = TRUE)
out_new <- data.frame(hex_id = hex_ids[fl[, "col"]], date = M$date[fl[, "row"]],
                      lst = X_raw[fl], deviation = D[fl], run_id = run_id)

n_valid <- sapply(files, \(f) global(!is.na(rast(f)), "sum")[1, 1])
inv_new <- data.frame(date = new_dates) |>
  left_join(data.frame(date = fdates, n_valid = as.integer(n_valid)), by = "date") |>
  mutate(n_valid   = coalesce(n_valid, 0L),
         pct_valid = round(100 * n_valid / ncell(ref), 1),
         status = case_when(!date %in% fdates ~ "missing composite",
                            n_valid == 0      ~ "empty (all cloud)",
                            pct_valid < 50    ~ "mostly cloud",
                            TRUE              ~ "usable"),
         run_id = run_id)

run_row <- data.frame(run_id = run_id, run_time = Sys.time(), run_type = "incremental update",
                      description = paste("Added MODIS", paste(yrs, collapse = ", "),
                                          "from staging; frozen v1.1 rules; earlier rows untouched"),
                      period_start = min(new_dates), period_end = max(new_dates),
                      method_version = "v1.1", n_lst_rows = nrow(lst_new))

# ---------- Section 7: write everything in ONE transaction (all or nothing) ----------
dbWithTransaction(con, {
  dbAppendTable(con, "run_log",       run_row)
  dbAppendTable(con, "lst_composite", inv_new)
  dbAppendTable(con, "lst_weekly",    lst_new)
  if (nrow(out_new)) dbAppendTable(con, "outlier_log", out_new)
})

# ---------- Section 8: what changed, and what did not ----------
cat("\nNew weeks by status:\n"); print(table(inv_new$status))
cat("\nNew hexagon-weeks by source:\n"); print(table(lst_new$source))

old_n_after <- q("SELECT COUNT(*) AS n, SUM(lst_filled) AS s FROM lst_weekly WHERE run_id = 1")
cat("\nEarlier rows (run 1) untouched:", isTRUE(all.equal(old_n_before, old_n_after)), "\n")

kids <- q("SELECT exam_date FROM child")
cat("Children whose exposure windows include the new weeks:",
    sum(kids$exam_date - 1 >= min(new_dates)), "(0 = exposure table needs no recalculation)\n")

print(sapply(dbListTables(con), \(t) q(paste("SELECT COUNT(*) AS n FROM", t))$n))
q("SELECT year(date) AS year, source, COUNT(*) AS n FROM lst_weekly GROUP BY ALL ORDER BY ALL") |> print()
q("SELECT run_id, run_type, period_start, period_end, n_lst_rows FROM run_log") |> print()
dbDisconnect(con)

# ---------- Section 9: tests AFTER the update ----------
cat("\n===== TESTS AFTER UPDATE =====\n")
test_file("tests/test_database.R")

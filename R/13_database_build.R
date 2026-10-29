# ============================================================
# 13_database_build.R - build the exposure database from scratch (run 1)
# ============================================================
library(DBI); library(duckdb); library(dplyr); library(sf)

dir.create("db", showWarnings = FALSE)
db_path <- "db/freiburg_exposure.duckdb"
if (file.exists(db_path)) file.remove(db_path)          # build = always from scratch
con <- dbConnect(duckdb(), dbdir = db_path)

# ---------- Section 1: schema (tables, types, keys) ----------
dbExecute(con, "CREATE TABLE run_log (
  run_id INTEGER PRIMARY KEY, run_time TIMESTAMP, run_type VARCHAR, description VARCHAR,
  period_start DATE, period_end DATE, method_version VARCHAR, n_lst_rows INTEGER)")

dbExecute(con, "CREATE TABLE lst_composite (
  date DATE PRIMARY KEY, n_valid INTEGER, pct_valid DOUBLE, status VARCHAR, run_id INTEGER)")

dbExecute(con, "CREATE TABLE hexagon (
  hex_id INTEGER PRIMARY KEY, inside_frac DOUBLE, green DOUBLE, cropland DOUBLE,
  built DOUBLE, blue DOUBLE, green_osm DOUBLE, blue_osm DOUBLE, ndvi DOUBLE,
  elev DOUBLE, geom_wkt VARCHAR)")

dbExecute(con, "CREATE TABLE lst_weekly (
  hex_id INTEGER, date DATE, lst_raw DOUBLE, lst_filled DOUBLE NOT NULL,
  source VARCHAR NOT NULL, method_version VARCHAR, run_id INTEGER,
  PRIMARY KEY (hex_id, date))")

dbExecute(con, "CREATE TABLE outlier_log (
  hex_id INTEGER, date DATE, lst DOUBLE, deviation DOUBLE, run_id INTEGER,
  PRIMARY KEY (hex_id, date))")

dbExecute(con, "CREATE TABLE child (
  child_id VARCHAR PRIMARY KEY, hex_id INTEGER, exam_date DATE, age DOUBLE, sex VARCHAR,
  ses_z DOUBLE, elev DOUBLE, green300 DOUBLE, ndvi300 DOUBLE, dist_park DOUBLE,
  park_access INTEGER)")

dbExecute(con, "CREATE TABLE exposure_lst (
  child_id VARCHAR, window_days INTEGER, coverage DOUBLE, lst_naive DOUBLE,
  lst_filled DOUBLE, run_id INTEGER, PRIMARY KEY (child_id, window_days))")

print(dbListTables(con))

# ---------- Section 2: load the data (run 1) ----------
hex <- st_read("data_cache/hex_green_blue.gpkg", quiet = TRUE)
hex_tab <- hex |>
  mutate(geom_wkt = st_as_text(st_geometry(hex))) |>         # shape stored as text (WKT)
  st_drop_geometry() |>
  select(any_of(c("hex_id", "inside_frac", "green", "cropland", "built", "blue",
                  "green_osm", "blue_osm", "ndvi", "elev", "geom_wkt")))

inv <- read.csv("data_cache/lst_inventory.csv") |>
  transmute(date = as.Date(date), n_valid = as.integer(n_valid), pct_valid, status, run_id = 1L)

lst_tab <- read.csv("data_cache/lst_hex_filled.csv") |>
  transmute(hex_id = as.integer(hex_id), date = as.Date(date), lst_raw = lst,
            lst_filled, source, method_version = "v1.1", run_id = 1L)

out_tab <- read.csv("outputs/tables/T1b_outliers_removed.csv") |>
  transmute(hex_id = as.integer(hex_id), date = as.Date(date), lst, deviation, run_id = 1L)

child_tab <- read.csv("data_cache/exposure_child_final.csv") |>
  transmute(child_id, hex_id = as.integer(hex_id), exam_date = as.Date(exam_date), age, sex,
            ses_z, elev, green300, ndvi300, dist_park, park_access = as.integer(park_access))

expo_tab <- read.csv("data_cache/exposure_final.csv") |>
  transmute(child_id, window_days = as.integer(window), coverage, lst_naive, lst_filled,
            run_id = 1L)

run1 <- data.frame(run_id = 1L, run_time = Sys.time(), run_type = "initial build",
                   description = "MODIS 2021-2023, gap-filling v1.1 (outlier screen + robust anomaly)",
                   period_start = min(inv$date), period_end = max(inv$date),
                   method_version = "v1.1", n_lst_rows = nrow(lst_tab))

dbAppendTable(con, "run_log", run1)
dbAppendTable(con, "lst_composite", inv)
dbAppendTable(con, "hexagon", hex_tab)
dbAppendTable(con, "lst_weekly", lst_tab)
dbAppendTable(con, "outlier_log", out_tab)
dbAppendTable(con, "child", child_tab)
dbAppendTable(con, "exposure_lst", expo_tab)

# ---------- Section 3: check counts, then ask questions in SQL ----------
counts <- sapply(dbListTables(con), \(t) dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM", t))$n)
print(counts)

# Provenance by year: how much of each year is observed vs filled?
dbGetQuery(con, "
  SELECT year(date) AS year, source, COUNT(*) AS n
  FROM lst_weekly GROUP BY ALL ORDER BY ALL") |> print()

# A join across three tables: annual heat of each child next to their hexagon's built-up share
dbGetQuery(con, "
  SELECT c.child_id, c.exam_date, h.built, ROUND(e.lst_filled, 2) AS lst_365
  FROM child c
  JOIN exposure_lst e USING (child_id)
  JOIN hexagon h USING (hex_id)
  WHERE e.window_days = 365
  ORDER BY h.built DESC LIMIT 5") |> print()

# ---------- Section 4: the rules work - try to break them ----------
dup <- dbGetQuery(con, "SELECT * FROM lst_weekly LIMIT 1")
r1 <- try(dbAppendTable(con, "lst_weekly", dup), silent = TRUE)
cat("Duplicate hexagon-week blocked:", inherits(r1, "try-error"), "\n")

bad <- dup; bad$lst_filled <- NA
bad$date <- as.Date("1999-01-01")
r2 <- try(dbAppendTable(con, "lst_weekly", bad), silent = TRUE)
cat("Missing filled value blocked:", inherits(r2, "try-error"), "\n")

dbDisconnect(con)
cat("Database size (MB):", round(file.size(db_path) / 1e6, 1), "\n")

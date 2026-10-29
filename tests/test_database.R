# ============================================================
# tests/test_database.R - automated integrity checks for the exposure database
# ============================================================
library(testthat); library(DBI); library(duckdb)
db <- file.path(rprojroot::find_root(rprojroot::has_file("renv.lock")),
                "db", "freiburg_exposure.duckdb")
con <- dbConnect(duckdb(), dbdir = db)
qry <- function(sql) dbGetQuery(con, sql)
n   <- function(sql) as.numeric(qry(sql)[1, 1])

test_that("one row per hexagon-week (no duplicates)", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly"),
               n("SELECT COUNT(*) FROM (SELECT DISTINCT hex_id, date FROM lst_weekly)"))
})

test_that("complete grid: every hexagon has every composite week", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly"),
               n("SELECT COUNT(*) FROM hexagon") * n("SELECT COUNT(*) FROM lst_composite"))
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly WHERE date NOT IN (SELECT date FROM lst_composite)"), 0)
})

test_that("no gaps left: every hexagon-week has a filled value", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly WHERE lst_filled IS NULL"), 0)
})

test_that("source labels are only the allowed ones", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly
                  WHERE source NOT IN ('observed', 'filled_clim_anomaly', 'removed_outlier')"), 0)
})

test_that("observed values are never altered by gap-filling", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly WHERE source = 'observed'
                  AND (lst_raw IS NULL OR abs(lst_filled - lst_raw) > 1e-9)"), 0)
})

test_that("every removed outlier is logged with its reason", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly w WHERE source = 'removed_outlier'
                  AND NOT EXISTS (SELECT 1 FROM outlier_log o
                                  WHERE o.hex_id = w.hex_id AND o.date = w.date)"), 0)
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly WHERE source = 'removed_outlier'"),
               n("SELECT COUNT(*) FROM outlier_log"))
})

test_that("every row traces back to a logged pipeline run", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly    WHERE run_id NOT IN (SELECT run_id FROM run_log)"), 0)
  expect_equal(n("SELECT COUNT(*) FROM lst_composite WHERE run_id NOT IN (SELECT run_id FROM run_log)"), 0)
})

test_that("temperatures are physically plausible (-40 to 70 °C)", {
  expect_equal(n("SELECT COUNT(*) FROM lst_weekly WHERE lst_filled < -40 OR lst_filled > 70"), 0)
})

test_that("land-cover shares stay between 0 and 1", {
  expect_equal(n("SELECT COUNT(*) FROM hexagon WHERE
                  green < 0 OR green > 1 OR built < 0 OR built > 1 OR cropland < 0 OR cropland > 1
                  OR blue < 0 OR blue > 1 OR inside_frac < 0 OR inside_frac > 1"), 0)
})

test_that("children: 500, unique, each living in a known hexagon", {
  expect_equal(n("SELECT COUNT(*) FROM child"), 500)
  expect_equal(n("SELECT COUNT(DISTINCT child_id) FROM child"), 500)
  expect_equal(n("SELECT COUNT(*) FROM child WHERE hex_id NOT IN (SELECT hex_id FROM hexagon)"), 0)
})

test_that("each child has exactly the 7-, 30- and 365-day windows, all filled", {
  expect_equal(n("SELECT COUNT(*) FROM (SELECT child_id FROM exposure_lst
                  GROUP BY child_id HAVING COUNT(*) <> 3)"), 0)
  expect_setequal(qry("SELECT DISTINCT window_days FROM exposure_lst")$window_days, c(7, 30, 365))
  expect_equal(n("SELECT COUNT(*) FROM exposure_lst WHERE lst_filled IS NULL"), 0)
})

test_that("composite calendar: each year starts on 1 January, steps of 8 days", {
  d  <- qry("SELECT date FROM lst_composite ORDER BY date")$date
  ok <- sapply(split(d, format(d, "%Y")),
               \(x) format(x[1], "%m-%d") == "01-01" && all(as.numeric(diff(x)) == 8))
  expect_true(all(ok))
})

dbDisconnect(con)

# Satellite exposure linkage for a child cohort — Freiburg, 2021–2024

**📄 [Read the full report](https://damanijay.github.io/Freiburg-exposure-linkage/report.html)**

A reproducible pipeline that turns raw satellite products into individual-level
environmental exposures, with validated gap-filling, provenance tracking and an
incremental update procedure.

## What it does

- Downloads and quality-screens MODIS land surface temperature (184 8-day composites, 2021–2024)
- Derives green space, cropland and built-up shares from ESA WorldCover, cross-checked against OpenStreetMap
- Computes 10 m NDVI from Sentinel-2 and elevation from Copernicus DEM
- Diagnoses missing data, then validates three gap-filling methods against withheld observations
- Links exposures to 500 simulated children over 7-, 30- and 365-day windows before each examination date
- Stores everything in DuckDB with full provenance, under 19 automated integrity tests

## Selected findings

| Finding | Result |
|---|---|
| Gaps are seasonal and clustered | 14.4% of hexagon-weeks missing; 46% in October, 0.9% in August |
| Ignoring gaps biases exposure | Annual exposure overestimated by 1.43 °C, same direction for nearly every child |
| Satellite tracks real weather | Deseasonalised r = 0.82 against DWD station Freiburg |
| Terrain confounds greenness | Apparent cooling per 10 pp green falls 45% after adjusting for elevation |
| Window choice matters | Exam month explains 92% of 7-day exposure; location explains 63% of 365-day exposure |

## Running it

```r
renv::restore()        # install pinned package versions
# then run R/ scripts in numerical order
```

Downloads are cached, so re-runs do not re-fetch. Satellite data and the database are
not in the repository; they are produced by the pipeline.

## Structure

```
R/                 numbered pipeline scripts
tests/             automated database integrity tests
outputs/tables/    result tables referenced in the report
report.Rmd         full methods and results
```

## Data sources

MODIS MOD11A2 v061 · Sentinel-2 L2A · ESA WorldCover 2021 · Copernicus DEM GLO-30 ·
OpenStreetMap · DWD Climate Data Center (station 01443)

## Author

Jay Damani — jaydamani210@gmail.com

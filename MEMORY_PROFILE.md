# Memory usage and profiling (Meta AML Shiny)

## Quick profile (simulated app load)

From the project root:

```bash
Rscript memory_profile_app.R
```

This loads the same data and libraries as the app (without starting Shiny) and reports memory (MB) at each stage.

## Results (typical run)

| Stage | Memory (MB) | Delta |
|-------|-------------|--------|
| 1. After gc(reset) | ~20 | baseline |
| 2. After shiny/ggplot2/survival/DT/scales | ~140 | +120 |
| 3. After optional packages (survminer, etc.) | ~261 | +120 |
| 4. After readRDS(AML_Meta_Cohort_v2.rds) | ~264 | +3 |
| 5–6. After normalize + temp CSV round-trip | ~264 | ~0 |
| 7. Filtered data (De novo, All) | ~267 | +2 |
| 8. Two cached filter results (analyses + meta_aml4) | ~267 | +0.6 |
| 9–10. Co-occurrence matrix + oncoprint prep | ~267 | negligible |
| 11. drug_sensitivity_precomputed.rds | ~269 | +2 |
| 12. load_beataml2() | ~288 | +19 |
| **Peak (max used)** | **~410** | |

## Where memory goes

- **Libraries (shiny, ggplot2, survival, DT, survminer, etc.)**: ~260 MB. Largest single contributor.
- **Cohort data**: Small (~3 MB for 13k rows × 23 cols). Grows with cohort size.
- **Filtered-data cache**: App keeps at most one cache per tab (Analyses and Meta AML4). When you switch main tab, the other tab’s cache is cleared and `gc()` is run so memory is returned sooner.
- **Beat AML2**: +~19 MB when drug data is loaded.
- **Precomputed drug**: +~2 MB if `drug_sensitivity_precomputed.rds` exists.

## Safeguards already in the app

1. **Single active cache**: Only the current main tab (Analyses or Meta AML4) holds a full filtered dataset; the other tab’s cache is cleared on switch.
2. **gc() on tab switch**: After clearing the other tab’s cache, `gc()` is called to return memory to the OS.
3. **gc() before Drug tab**: When opening the Single Gene “Drug Sensitivity” sub-tab, `gc()` runs before loading drug content to reduce peak memory.
4. **Precomputed drug**: Using `drug_sensitivity_precomputed.rds` avoids recomputing drug sensitivity at runtime.

## Recommendations

1. **Deploy / resource limits**: A 512 MB–1 GB RAM limit per process is typically enough for normal use (libraries + data + one cached filter + Beat AML + overhead). Peak in testing was ~410 MB.
2. **Large cohorts**: If the cohort RDS grows (e.g. 100k+ rows), consider limiting genes in co-occurrence and oncoprint (e.g. top N by frequency) to keep matrix sizes bounded.
3. **Temp CSV**: The app currently writes the normalized cohort to a temp CSV and reads it back. This was kept for compatibility; removing it would avoid a brief duplication during startup.
4. **Re-profile after changes**: Run `Rscript memory_profile_app.R` after adding large data or new packages to confirm memory stays acceptable.

## Optional: profile live app

To see process memory while using the app (e.g. on macOS):

```bash
# Start app in background, then in another terminal:
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)' &
sleep 10
ps -o rss= -p $(pgrep -f "runApp") | awk '{print $1/1024 " MB"}'
```

Or use the system monitor (Activity Monitor / Task Manager) and filter by R or the browser process.

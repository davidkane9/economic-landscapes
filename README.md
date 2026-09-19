# Economic Landscapes

Code and results for the paper

> David Kane, "Economic Landscapes: Tunable Complexity and the Limits of
> Organizational Search."

A firm allocates a fixed integer budget among *N* inputs. Performance is the
quadratic form *x'Qx* over the budget simplex, and the firm hill-climbs among
neighboring allocations defined by a one-dollar transfer operator. The paper
maps how the difficulty of that search varies with the curvature of *Q* and
with the connection density *C/N* of the search neighborhood.

Everything needed to reproduce the simulations, figures, and tables is here.
The study uses no external data.

## Contents

| Path | What it is |
|------|------------|
| `src/landscape.cpp` | The simulation engine (C++17, no dependencies). Its header comment documents the model and every command. |
| `sim/run_sweeps.sh` | Runs every simulation in the paper, with fixed seeds, and writes the CSVs in `sim/`. |
| `sim/*.csv` | The simulation results the paper's figures and tables are built from. |
| `R/simulations.R` | R reference implementation of the model, the stylized landscapes of Fig. 1, and the local-optima count. |
| `R/complexity.R` | Figures 2–6 and Table 1, drawn from `sim/*.csv`. |
| `R/eigen_fractions.R` | The negative-eigenvalue fractions quoted in the text and in the caption of Fig. 2. |
| `economic-landscapes.qmd` | The manuscript (Quarto). |
| `economic-landscapes.pdf` | The rendered manuscript. |
| `author.yaml` | Author metadata. Without it the manuscript renders in the anonymized form used for double-anonymous review. |
| `references.bib`, `springer-basic-author-date.csl` | Bibliography and citation style. |

An ODD-protocol description of the model is the appendix of the manuscript.

## Reproducing the paper

Requirements: a C++17 compiler; R (4.x) with `ggplot2`, `dplyr`, `tidyr`, and
`knitr`; [Quarto](https://quarto.org) (1.4 or later) with a LaTeX
distribution for the PDF (`quarto install tinytex`).

1. **Simulations** (optional; the results are already in `sim/`):

   ```
   bash sim/run_sweeps.sh
   ```

   This compiles the engine and runs about eighty-five jobs in parallel (see
   the note on run time below). Set `CXX` to choose a compiler.

2. **Eigenvalue fractions** quoted in the text:

   ```
   Rscript R/eigen_fractions.R
   ```

3. **Manuscript**, with all figures and tables:

   ```
   quarto render economic-landscapes.qmd --metadata-file author.yaml
   ```

   The figures and tables read `sim/*.csv`; only the stylized landscapes and
   the local-optima sample are computed at render time (two to three minutes).

## How closely a re-run reproduces the stored results

The seeds are fixed in `sim/run_sweeps.sh`, and on a given toolchain the output
is deterministic. Re-running jobs on the toolchain below gave:

| File | Re-run vs stored |
|------|------------------|
| `transition`, `regime`, `budget`, `operator`, `enumerate`, `qplane`, `lon`, `autocorr` | identical, or within 0.0001 (the last printed digit; floating-point summation order) |
| `plane`, `diag` | agree within sampling error (differences of about one standard error). These two files were written by an earlier build of the engine that consumed the random stream differently, so the same seeds now give different landscapes. The conclusions drawn from them do not change. |

The engine draws random numbers with the C++ standard library's
`uniform_real_distribution`, whose output is not specified across
standard-library implementations. On a different toolchain expect agreement
within the reported standard errors, not digit for digit.

**Run time.** Most jobs finish in seconds to minutes. The (N, C) plane does
not: its jobs for N above 30 each run for many hours, so the full script takes
the better part of a day on a ten-core laptop. To check the pipeline quickly,
run single jobs by hand; the engine's header comment lists every command, for
example

```
c++ -O3 -std=c++17 -o landscape src/landscape.cpp
./landscape sweep 10 50 uniform 1,2,3,4,5,6,7,8,9 1000 8 101   # ~5 s; rows of sim/transition.csv
```

Stored results: Apple clang 21 on macOS (arm64), R 4.5, Quarto 1.9.

## License

MIT. See `LICENSE`.

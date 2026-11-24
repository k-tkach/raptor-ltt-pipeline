# raptor-ltt-pipeline

Reproducible R pipeline to build **lineage-through-time (LTT) curves** and detect
**shifts in the rate of lineage accumulation** across a posterior sample of
phylogenies, accounting for phylogenetic uncertainty.

It was developed for my PhD research on the diversification of the four raptor
orders (Accipitriformes, Falconiformes, Strigiformes and Cathartiformes) using
200 trees from [BirdTree.org](https://birdtree.org/), but it works with any clade
and any set of ultrametric trees.

## What it does

1. Prunes every tree to the species of interest.
2. Builds the **maximum clade credibility (MCC) tree** from the posterior sample.
3. Computes **LTT curves** for each clade in every tree, on a common time grid.
4. Summarises **crown ages** (median and 95% interval across trees).
5. Fits **piecewise (segmented) regressions** of log(lineages) vs. time in every
   tree, with 0–2 breakpoints selected by **BIC**, and reports the age of each
   breakpoint and the slopes before and after it.
6. Summarises how consistently each breakpoint appears across trees, so results
   do not depend on a single phylogeny.

## Repository structure

```
R/ltt_functions.R       core functions (pruning, LTT, breakpoint detection)
R/run_ltt_analysis.R    full analysis on real data (BirdTree + taxonomy)
R/example_simulated.R   self-contained example on simulated trees (~3 s)
data/README.md          how to obtain the input data
```

## Quick start

Requires R (≥ 4.1) and the packages `ape`, `dplyr`, `ggplot2` and, optionally,
`phangorn` (for the MCC tree).

```r
install.packages(c("ape", "dplyr", "ggplot2", "phangorn"))
```

Test the pipeline on simulated data (no downloads needed), from the repository root:

```bash
Rscript R/example_simulated.R
```

Run the full analysis after placing the input files in `data/` (see `data/README.md`):

```bash
Rscript R/run_ltt_analysis.R
```

## Outputs (`output/`)

| File | Content |
|---|---|
| `crown_ages.csv` | Crown age per clade: median, 95% interval and MCC value |
| `breakpoints_per_tree.csv` | Breakpoint ages and slopes for every tree and the MCC tree |
| `breakpoint_frequency.csv` | Proportion of trees with a breakpoint in each 5-Ma bin |
| `ltt_by_clade.png` | LTT curves of all trees (grey) and the MCC tree (black) |
| `breakpoints_histogram.png` | Distribution of breakpoint ages across trees |
| `sessionInfo.txt` | R and package versions, for reproducibility |

## Methodological notes

- Breakpoints are fitted on log(N) to compare per-lineage rates between periods.
- The analysis is **exploratory**: successive nodes are not independent, so BIC
  values are used to compare models, not as formal hypothesis tests.
- Recent slowdowns in LTT curves can reflect incomplete sampling or taxonomy
  ("pull of the present"); interpret breakpoints near the present with caution.

## Reference

Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012).
The global diversity of birds in space and time. *Nature*, 491, 444–448.

## Author

**Kevin A. Tkach** – PhD candidate in Ecology, University of Buenos Aires
[LinkedIn](https://www.linkedin.com/in/kevintkach)

## License

MIT

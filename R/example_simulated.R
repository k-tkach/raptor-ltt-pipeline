# =============================================================================
# example_simulated.R
# Self-contained example: runs the full pipeline on simulated trees so the
# code can be tested in under a minute without downloading any data.
#
# Run from the repository root:  Rscript R/example_simulated.R
# =============================================================================

source("R/ltt_functions.R")
set.seed(42)

# 50 simulated ultrametric trees with 120 tips scaled to a 70-Ma root
n_tips <- 120
tips   <- paste0("Sp_", seq_len(n_tips))
trees  <- lapply(1:50, function(i) {
  tr <- rcoal(n_tips, tip.label = sample(tips))
  tr$edge.length <- tr$edge.length * 70 / max(branching.times(tr))
  tr
})
class(trees) <- "multiPhylo"

# Three arbitrary "clades" of different size
taxonomy <- data.frame(species = tips,
                       clade = rep(c("Clade_A", "Clade_B", "Clade_C"),
                                   times = c(60, 40, 20)))

res <- run_ltt_pipeline(trees, taxonomy, age_grid = seq(0, 75, by = 0.5), kmax = 2)

print(res$crown_ages)
print(head(res$breakpoints))
stopifnot(nrow(res$crown_ages) == 3,
          all(c("tree", "k", "age_Ma") %in% names(res$breakpoints)))
message("Example finished successfully.")

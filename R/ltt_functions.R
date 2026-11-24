# =============================================================================
# ltt_functions.R
# Functions to build lineage-through-time (LTT) curves across a posterior
# sample of phylogenies and to detect shifts in the rate of lineage
# accumulation with piecewise (segmented) regression.
#
# Author: Kevin A. Tkach (University of Buenos Aires)
# License: MIT
# =============================================================================

library(ape)
library(dplyr)

#' Prune every tree in a multiPhylo object to a set of species.
prune_trees <- function(trees, species) {
  pruned <- lapply(trees, function(tr) keep.tip(tr, intersect(species, tr$tip.label)))
  class(pruned) <- "multiPhylo"
  pruned
}

#' Node (branching) ages of one clade in one tree, oldest first.
#' Returns NULL when the clade has fewer than 3 species in the tree.
clade_node_ages <- function(tree, species) {
  spp <- intersect(species, tree$tip.label)
  if (length(spp) < 3) return(NULL)
  sort(branching.times(keep.tip(tree, spp)), decreasing = TRUE)
}

#' Number of lineages alive at each age of a time grid.
#' Ages older than the crown age are returned as NA.
lineages_on_grid <- function(node_ages, age_grid) {
  n <- 1 + length(node_ages) - findInterval(age_grid, sort(node_ages))
  n[age_grid >= max(node_ages)] <- NA
  n
}

#' Piecewise linear regression of log(N lineages) on time, with 0 to kmax
#' breakpoints selected by BIC. The fit is continuous at the breakpoints.
#'
#' Exploratory: successive nodes are not independent observations, so BIC
#' values should not be read as formal hypothesis tests.
fit_breakpoints <- function(time, n, kmax = 2, min_pts = 5, max_cand = 40) {
  x <- time; y <- log(n); N <- length(x)
  if (N < 2 * min_pts + 1) kmax <- 0
  idx  <- seq(min_pts, max(min_pts, N - min_pts))
  cand <- unique(x[idx[round(seq(1, length(idx),
                                 length.out = min(max_cand, length(idx))))]])

  design <- function(b) {
    X <- cbind(1, x)
    if (length(b)) X <- cbind(X, vapply(b, function(bi) pmax(x - bi, 0), numeric(N)))
    X
  }
  fit_one <- function(b) {
    fit <- lm.fit(design(b), y)
    rss <- sum(fit$residuals^2)
    list(b = b, coef = fit$coefficients,
         bic = N * log(rss / N) + (2 + 2 * length(b)) * log(N))
  }

  models <- list(fit_one(numeric(0)))
  if (kmax >= 1) models <- c(models, lapply(cand, fit_one))
  if (kmax >= 2 && length(cand) >= 2) {
    pairs <- t(combn(cand, 2))
    pairs <- pairs[(pairs[, 2] - pairs[, 1]) > 0.1 * diff(range(x)), , drop = FALSE]
    models <- c(models, lapply(seq_len(nrow(pairs)), function(i) fit_one(pairs[i, ])))
  }
  models[[which.min(vapply(models, `[[`, numeric(1), "bic"))]]
}

#' Breakpoints (age in Ma) and slopes before/after, for one set of node ages.
breakpoints_table <- function(node_ages, tree_id, kmax = 2) {
  node_ages <- sort(node_ages, decreasing = TRUE)
  m <- fit_breakpoints(-node_ages, seq_along(node_ages) + 1, kmax = kmax)
  if (!length(m$b)) {
    return(data.frame(tree = as.character(tree_id), k = 0, age_Ma = NA,
                      slope_before = NA, slope_after = NA))
  }
  slopes <- cumsum(m$coef[-1])            # slope of each segment
  data.frame(tree = as.character(tree_id), k = length(m$b), age_Ma = unname(-m$b),
             slope_before = unname(slopes[seq_along(m$b)]),
             slope_after  = unname(slopes[seq_along(m$b) + 1]),
             row.names = NULL)
}

#' Full analysis for several clades across a posterior sample of trees.
#'
#' @param trees     multiPhylo object (e.g., a posterior sample from BirdTree).
#' @param taxonomy  data.frame with columns `species` (Genus_species) and `clade`.
#' @param age_grid  ages (Ma before present) at which LTT curves are evaluated.
#' @param kmax      maximum number of breakpoints per tree.
#' @return list with LTT curves, crown ages, breakpoints and the MCC tree.
run_ltt_pipeline <- function(trees, taxonomy, age_grid = seq(0, 80, by = 0.25),
                             kmax = 2) {
  stopifnot(inherits(trees, "multiPhylo"),
            all(c("species", "clade") %in% names(taxonomy)))
  clades <- sort(unique(taxonomy$clade))
  spp_by_clade <- split(taxonomy$species, taxonomy$clade)

  trees <- prune_trees(trees, taxonomy$species)
  message("Trees: ", length(trees), " | tips per tree: ", Ntip(trees[[1]]))

  # Maximum clade credibility tree (requires phangorn)
  mcc <- if (requireNamespace("phangorn", quietly = TRUE)) {
    phangorn::maxCladeCred(trees)
  } else {
    warning("phangorn not installed: tree 1 is used instead of the MCC tree")
    trees[[1]]
  }

  ages     <- setNames(lapply(clades, function(cl)
                lapply(trees, clade_node_ages, species = spp_by_clade[[cl]])), clades)
  ages_mcc <- setNames(lapply(clades, function(cl)
                clade_node_ages(mcc, spp_by_clade[[cl]])), clades)

  curves <- bind_rows(lapply(clades, function(cl) bind_rows(lapply(
    seq_along(ages[[cl]]), function(i)
      data.frame(tree = i, clade = cl, age = age_grid,
                 n = lineages_on_grid(ages[[cl]][[i]], age_grid))))))
  curve_mcc <- bind_rows(lapply(clades, function(cl)
    data.frame(clade = cl, age = age_grid, n = lineages_on_grid(ages_mcc[[cl]], age_grid))))

  crown <- bind_rows(lapply(clades, function(cl)
    data.frame(clade = cl, crown_age = vapply(ages[[cl]], max, numeric(1))))) %>%
    group_by(clade) %>%
    summarise(median = median(crown_age), q025 = quantile(crown_age, 0.025),
              q975 = quantile(crown_age, 0.975), .groups = "drop") %>%
    mutate(mcc = vapply(clade, function(cl) max(ages_mcc[[cl]]), numeric(1)))

  breaks <- bind_rows(lapply(clades, function(cl) {
    bind_rows(c(Map(breakpoints_table, ages[[cl]], seq_along(ages[[cl]]), kmax = kmax),
                list(breakpoints_table(ages_mcc[[cl]], "MCC", kmax = kmax)))) %>%
      mutate(clade = cl, .before = 1)
  }))

  list(curves = curves, curve_mcc = curve_mcc, crown_ages = crown,
       breakpoints = breaks, mcc_tree = mcc)
}

# =============================================================================
# run_ltt_analysis.R
# Lineage-through-time (LTT) curves and diversification-rate breakpoints for
# raptor orders (Accipitriformes, Falconiformes, Strigiformes, Cathartiformes)
# across a posterior sample of bird phylogenies.
#
# Inputs (see data/README.md):
#   data/birdtree_sample.tre  - multiPhylo, e.g. 200 trees from birdtree.org
#   data/taxonomy.csv         - columns: species (Genus_species), clade
#
# Run from the repository root:  Rscript R/run_ltt_analysis.R
# =============================================================================

source("R/ltt_functions.R")
library(ggplot2)

# ---- Settings -----------------------------------------------------------
tree_file     <- "data/birdtree_sample.tre"
taxonomy_file <- "data/taxonomy.csv"
age_grid      <- seq(0, 80, by = 0.25)   # Ma before present
kmax          <- 2                       # maximum breakpoints per tree
out_dir       <- "output"
dir.create(out_dir, showWarnings = FALSE)
set.seed(1)

# ---- Data ---------------------------------------------------------------
trees    <- read.tree(tree_file)
taxonomy <- read.csv(taxonomy_file, stringsAsFactors = FALSE)

# ---- Analysis -----------------------------------------------------------
res <- run_ltt_pipeline(trees, taxonomy, age_grid = age_grid, kmax = kmax)

write.csv(res$crown_ages,  file.path(out_dir, "crown_ages.csv"),  row.names = FALSE)
write.csv(res$breakpoints, file.path(out_dir, "breakpoints_per_tree.csv"), row.names = FALSE)
print(res$crown_ages)

# Share of trees with a breakpoint in each 5-Ma bin (robustness across trees)
summary_breaks <- res$breakpoints %>%
  filter(tree != "MCC", !is.na(age_Ma)) %>%
  mutate(bin_Ma = floor(age_Ma / 5) * 5) %>%
  group_by(clade, bin_Ma) %>%
  summarise(n_trees = n_distinct(tree), .groups = "drop") %>%
  mutate(prop_trees = n_trees / length(trees))
write.csv(summary_breaks, file.path(out_dir, "breakpoint_frequency.csv"), row.names = FALSE)

# ---- Figures ------------------------------------------------------------
p_ltt <- ggplot() +
  geom_line(data = res$curves, aes(-age, n, group = tree),
            colour = "grey60", alpha = 0.2, na.rm = TRUE) +
  geom_line(data = res$curve_mcc, aes(-age, n),
            colour = "black", linewidth = 0.8, na.rm = TRUE) +
  scale_y_log10() +
  facet_wrap(~clade, scales = "free_y") +
  labs(x = "Time (Ma)", y = "Number of lineages (log scale)",
       caption = "Grey: posterior trees; black: maximum clade credibility tree") +
  theme_classic()
ggsave(file.path(out_dir, "ltt_by_clade.png"), p_ltt, width = 10, height = 7, dpi = 300)

p_breaks <- ggplot(filter(res$breakpoints, tree != "MCC", !is.na(age_Ma)), aes(age_Ma)) +
  geom_histogram(binwidth = 2, fill = "steelblue", colour = "white") +
  scale_x_reverse() +
  facet_wrap(~clade, scales = "free_y") +
  labs(x = "Breakpoint age (Ma)", y = "Number of trees") +
  theme_classic()
ggsave(file.path(out_dir, "breakpoints_histogram.png"), p_breaks, width = 10, height = 7, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Results in ./", out_dir)

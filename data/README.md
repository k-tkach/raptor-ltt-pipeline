# Input data

The raw inputs are **not** included in this repository because they are
distributed by their original providers. To reproduce the analysis:

1. **Phylogenies** – download a posterior sample of bird trees (e.g., 200 trees,
   Hackett backbone) from [BirdTree.org](https://birdtree.org/) (Jetz et al. 2012)
   and save it as `data/birdtree_sample.tre` (Newick, multiPhylo).
2. **Taxonomy** – create `data/taxonomy.csv` with two columns:

   | species | clade |
   |---|---|
   | Falco_sparverius | Falconiformes |
   | Buteo_magnirostris | Accipitriformes |

   Species names must match the tip labels of the trees (`Genus_species`).
   Any taxonomy can be used (e.g., the HBW and BirdLife International checklist).

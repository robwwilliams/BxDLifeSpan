## BxDLifeSpan

## Rob Williams's fork (21 November 2026)

Danny Arends's code, forked for the Nature Aging paper. Where the mapping stands, why each model was chosen, and the results so far: [BXD_mapping_status_2026-11-23.md](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/notes/BXD_mapping_status_2026-11-23.md) in bxd_mortality_rates_2026. Danny's files are unchanged; two files are added.

- **`actuary_map.R`** draws Danny's actuarial heat maps with one fixed geometry, so that panels from every script superimpose with pixel accuracy.
  - It keeps his colors, bins, row lines and chromosome boxes.
  - **Margins are set in inches**, the same for every panel. The right margin always has room for 4-digit Ns. Before, `QTL_plots.R` used R's default right margin (2.1 lines) and `interaction_plot.R` used 4.1 lines, so the chromosome area of the strain-means maps was 0.4 inch wider than that of the individual maps.
  - **The rows are always `plotTPS`** (50 to 770 days), and the chromosome layout always comes from the full map, with 15 Mb gaps.
  - **Every page has a white background** covering it, so the bounding box Illustrator sees is the page itself.
  - **Runs of neighboring markers with the same color** are drawn as one rectangle: the same image, without hairline seams.
  - **Output is SVG** (svglite, Arial, each label one text element), or PDF if the file name ends in `.pdf`.
- **`figure_maps.R`** draws every map of Fig. 2 with it, in `output/maps/`:
  - one SVG per panel, 24 × 5.5 inches;
  - Danny's pages of two;
  - all panels on one page.

  It reads Danny's result files in `output/` when they are present. Otherwise it reads the two workbooks he sent on 5 November 2026:

  ```bash
  Rscript figure_maps.R QTL_Lifespan_05Nov2026.xlsx QTL_BW_5Nov.xlsx
  ```

  The numbers on the right axis are counted with Danny's rules. They match his figure: 93 to 46 strains on chow, 87 to 14 on high fat, and 1,745 to 402 mice.

  Each panel prints its model, the mapping script and its commit, and the results file, in small type under the map. Titles are left-justified and in plain words. The legend sits in the right margin, clear of the maps, and says which allele each color favors: warm colors (yellow to red) for *D*, cool colors (pale green to dark blue) for *B*. Ticks mark every 25 Mb, longer every 50 Mb. Type sizes follow Rob's working Illustrator file.
- **`Figure_legends.md`** holds the Nature-style legends drafted for Figures 2 and 3, and **`Methods_draft.md`** the matching Methods paragraphs.
- **`ctl_signs.R`** reports the sign of the CTL correlations by diet and age, for the Figure 3 legend.
- **`compare_kinship.R`** sets Danny's individual-mouse maps beside the R/qtl2 kinship maps, T-age by T-age.
- **`qtl2_strainmeans.R`** maps strain means with R/qtl2 and LOCO kinship under three weightings: equal (Rob's default), mice per strain (Danny's) and inverse variance. **`compare_weights.R`** compares them. **`qtl2_cross.R`** builds the cross and kinship for both R/qtl2 scripts.
- **`qtl2_strainF1.R`** maps BXD strain means and BXD × BXD F1 means together (equal weights, at least 3 mice), with F1 type as a cofactor, for additive and dominance effects. F1 genotype probabilities come from both parents. It also lists where each F1 type is heterozygous on proximal Chr 1 (`output/strainF1_chr1.txt`).
- **`chr1_zoom.R`** zooms in on proximal Chr 1 (0 to 50 Mb), strain means alone against strain and F1 means together: zoomed actuarial maps (`output/maps/Chr1_zoom_maps.svg`), profiles at the peak T-ages over where each F1 type is heterozygous (`output/maps/Chr1_zoom_profiles.svg`), and the 1.5-LOD support interval at every T-age (`output/chr1_zoom_intervals.txt`). **`chr1_f1_dropone.R`** repeats the strain + F1 scan of Chr 1 with each F1 type, and each class of F1, left out (`output/chr1_f1_dropone.txt`).
- **`qtl2_strainF1_perm.R`** sets genome-wide thresholds for the strain + F1 maps (additive and dominance) from 1,000 permutations: strains among strains, F1 types within class. **`significant_loci_F1.R`** lists the loci that pass (`output/significant_loci_F1.txt`). **`chr1_twolocus.R`** asks whether high-fat proximal Chr 1 holds one locus or two, with a two-locus scan and a parametric bootstrap (`output/chr1_twolocus.txt`, `output/maps/Chr1_twolocus.svg`). **`strainF1_setup.R`** builds the strain + F1 units, genotypes and kinship for these scripts, as in `qtl2_strainF1.R`.
- **`chr1_recombinants.R`** lists the BXD strains with B/D switches in the chow block (4.86-5.32 Mb) and the high-fat interval (34.4-38.7 Mb), with their numbers of mice with lifespans (`output/chr1_recombinants.txt`): the strains, and F1s, that could split the intervals.
- **`qtl2_individualsF1.R`** maps individual BXD and F1 mice together (chow, high fat, both diets, genotype × diet), with F1 class and, for both diets, F1 class × diet as cofactors; Rob's Fig. 2 of 23 November 2026 (`output/maps/Figure2_lifespan_F1.svg`, drawn by `figure_maps.R`). Its thresholds come from `QTL/qtl2_individualsF1_perm.R` in bxd_mortality_rates_2026, run on ATHENA.
- **`cauchy_strainmeans.R`** combines each marker's P values over T-ages (the Cauchy combination of Arends et al. 2026) for the strain-mean maps with and without F1s, with genome-wide thresholds from 1,000 permutations (`output/cauchy_strainmeans.txt`).
- **`umhet3_maps.R`** draws the UM-HET3 actuarial maps of Arends et al. 2026 in this fork's exact geometry: a maternal map (*B* − *C*) and a paternal map (*D* − *H*), each on 1 d.f., by sex, with thresholds from 1,000 permutations within site and the Cauchy maxima (`output/umhet3_thresholds_<sex>.txt`). The BXD colors keep their meaning (warm *D*, cool *B*), *C* and *H* are gray, and the scale runs to >8. Pages: T50–T770 (the BXD rows), T50–T1100, and rows matched one to one with the BXD rows, on a 40-day grid (T45 to T1005) or by fraction alive (`output/maps/UMHET3_*`). It reads Supplementary Tables 7 and 9 from `../bxd_mortality_rates_2026/data/local/`.
- **`f1_heterosis.R`** tabulates heterosis (F1 mean minus midparent) for each F1 type, and tests whether it differs within a class by more than sampling error (`output/f1_heterosis.txt`).
- **`qtl2_permutations.R`** sets genome-wide thresholds for the strain-means maps by permuting strains (each strain takes another's whole actuarial series). **`significant_loci.R`** lists the loci that pass.
- **`qtl2_individuals_perm.R`** sets genome-wide thresholds for the individual-mouse maps by permuting whole strains within two strata (BXD1–42 and BXD43 on), with each strain's mice keeping their covariates. One permutation takes about 29 minutes on one worker. It saves each permutation as it finishes and resumes where it stopped. To run it on another machine (R with qtl2 and readxl):

  ```bash
  git clone https://github.com/robwwilliams/BxDLifeSpan.git
  git clone https://github.com/robwwilliams/bxd_mortality_rates_2026.git
  cd BxDLifeSpan && Rscript qtl2_individuals_perm.R 200 14 ../bxd_mortality_rates_2026/data/original/BDL_10002.xlsx
  ```

  The output, `output/qtl2_individuals_perm_maxima.txt`, comes back to this folder.
- **`run_CTL.R`** runs the first part of `bodyweight.R` (the interpolated weights), then `CTL.R`, so that the CTL panels can be drawn.
- **`qtl2_individuals.R`** maps individual mice with R/qtl2: linear mixed models with LOCO kinship, at every T-age.
  - It uses the same mice, genotypes and cofactors as `interactions.R`, with the kinship matrix in place of the random strain intercept.
  - Maps: chow alone, high fat alone, both diets (main effect), and genotype × diet.
  - `figure_maps.R` draws them as extra panels, with `Figure2_lifespan_kinship.svg`.
  - Rob: "even better is to simply map the individuals with a full kinship correction."

**Results on GitHub.** Danny's `.gitignore` keeps `output/` out of git. The small result tables and the main figures are added to git one by one (`git add -f`), so that the status note can link to them; the large result matrices (about 11 MB each) stay local.

Run order, from the fork's root:

```bash
Rscript run_CTL.R
Rscript qtl2_individuals.R
Rscript figure_maps.R QTL_Lifespan_05Nov2026.xlsx QTL_BW_5Nov.xlsx
```

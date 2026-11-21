## BxDLifeSpan

## Rob Williams's fork (21 November 2026)

Danny Arends's code, forked for the Nature Aging paper. Danny's files are unchanged; two files are added.

- **`actuary_map.R`** draws Danny's actuarial heat maps with one fixed geometry, so that panels from every script superimpose with pixel accuracy.
  - It keeps his colours, bins, row lines and chromosome boxes.
  - **Margins are set in inches**, the same for every panel. The right margin always has room for 4-digit Ns. Before, `QTL_plots.R` used R's default right margin (2.1 lines) and `interaction_plot.R` used 4.1 lines, so the chromosome area of the strain-means maps was 0.4 inch wider than that of the individual maps.
  - **The rows are always `plotTPS`** (50 to 770 days), and the chromosome layout always comes from the full map, with 15 Mb gaps.
  - **Every page has a white background** covering it, so the bounding box Illustrator sees is the page itself.
  - **Runs of neighbouring markers with the same colour** are drawn as one rectangle: the same image, without hairline seams.
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

  Each panel prints its model, the mapping script and its commit, and the results file, in small type under the map. Titles are left-justified and in plain words. The legend sits in the right margin, clear of the maps, and says what each colour means.
- **`run_CTL.R`** runs the first part of `bodyweight.R` (the interpolated weights), then `CTL.R`, so that the CTL panels can be drawn.
- **`qtl2_individuals.R`** maps individual mice with R/qtl2: linear mixed models with LOCO kinship, at every T-age.
  - It uses the same mice, genotypes and cofactors as `interactions.R`, with the kinship matrix in place of the random strain intercept.
  - Maps: chow alone, high fat alone, both diets (main effect), and genotype × diet.
  - `figure_maps.R` draws them as extra panels, with `Figure2_lifespan_kinship.svg`.
  - Rob: "even better is to simply map the individuals with a full kinship correction."

Run order, from the fork's root:

```bash
Rscript run_CTL.R
Rscript qtl2_individuals.R
Rscript figure_maps.R QTL_Lifespan_05Nov2026.xlsx QTL_BW_5Nov.xlsx
```

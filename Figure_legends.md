# Figure legends: draft for the Nature Aging paper

*Draft by Claude (Opus 5.5) for Rob Williams, 21 November 2026; Fig. 2 with the F1s added 23 November 2026. Not yet reviewed. The style follows Arends et al. 2026 (*Nature* 654:437), "Fig. N |" and bold panel letters, with American spelling, as *Nature Aging* uses (Rob, 24 November 2026). Each number comes from the files named under its figure.*

*Before submission:*
- *A Nature legend title states the finding. These titles describe instead. Genome-wide thresholds now exist for the strain-mean maps with and without F1s (1,000 permutations each) and for the individual-mouse maps without F1s (200 permutations, ATHENA); those with F1 mice are queued on ATHENA, and the body-weight and CTL maps still need theirs.*
- *The small model notes under each panel belong in these legends and in Methods; drop them from the final panels.*

## Fig. 2: lifespan, with the F1 hybrids (Rob's layout of 23 November 2026)

*Panels: `output/maps/Figure2_lifespan_F1.svg`. Rob: "a. strain means with F1s separately for both diets (first dyad); b. individuals with F1s stratified by CD and HFD (second dyad); c. individuals with F1s, all combined, all cofactors, main effect only; d. individuals all combined, diet x locus interaction effect only". His dyads a and b are panels **a**,**b** and **c**,**d** here; his c and d are **e** and **f**. The F1 mice are those of the table of 24 November 2026 (319 natural deaths). Thresholds for **c**–**f** wait for permutations with these F1 mice; the jobs now queued on ATHENA use the December 2024 table (308 F1 mice). The version below without F1s becomes Extended Data.*

**Fig. 2 | Actuarial genetic maps of lifespan in female BXD mice and BXD × BXD F1 hybrids on chow and high-fat diets.** Each map is a stack of genome scans, one per survivorship: the mice alive at a minimum age (T-age, left *y* axis), from T50 to T770 in 30-day steps (computed from T20), labeled every 60 days. Columns are the 21,056 markers of the BXD.8 genotypes on chromosomes 1–19 and X, by physical position. Each scan is a linear mixed model (R/qtl2) with a polygenic effect whose kinship matrix, computed over strains and F1 types together, leaves out the chromosome scanned (LOCO). Each of the 13 F1 types (6 short × short and 7 long × long crosses) takes its genotype probabilities from both parents, and genotype enters as the additive dosage of the *D* allele. F1 class (inbred, short × short, long × long) is a covariate throughout, so that the average hybrid vigor of a class is not read as a genotype effect.

Shading gives −log10 *P* of the genotype term (likelihood-ratio test, 1 d.f.), in bins of 1–2, 2–3, 3–4, 4–5 and >5; white is below 1. Warm colors mark markers at which the *D* allele is associated with longer life, and cool colors the *B* allele. Ticks under each chromosome mark every 25 Mb, longer every 50 Mb. Numbers of units (**a**,**b**) or mice (**c**–**f**) in each survivorship are on the right *y* axis.

**a**,**b**, Mean age at death of each BXD strain and F1 type with at least 3 mice alive at the T-age, on chow (**a**; 93 strains and 13 F1 types at T50, 46 and 13 at T770) and high fat (**b**; 87 and 13 at T50, 14 and 6 at T770), every unit weighted equally.
**c**,**d**, Age at death of individual mice on chow (**c**; 1,085 mice, 157 of them F1, at T50; 391 and 99 at T770) and high fat (**d**; 979 and 162 at T50; 155 and 45 at T770), with vivarium and set-up age (180-day classes) as further covariates.
**e**, Both diets, main effect of genotype, with diet and F1 class × diet as further covariates (2,064 mice, 319 F1, at T50; 546 and 144 at T770).
**f**, Genotype × diet interaction in the mice of **e**: warm colors where the *D* − *B* difference in lifespan is larger on high fat. F1 class × diet is a covariate so that the difference in hybrid vigor between the diets is not read as genotype × diet.

**Genome-wide significance.** In **a**,**b**, thresholds come from 1,000 permutations in which each strain took the whole actuarial series of another strain and each F1 type that of another F1 type of its class. At 5%, −log10 *P* for a single T-age is 4.7 on chow (T20–T740) and 4.6 on high fat (T20–T590); for any T-age, 5.8 and 5.5. Passing at a single T-age: on chow, Chr 1 near 5 Mb (*D* allele, 6.71 at T770; T770 lies one T-age past the stable range of the permutations, but 6.71 exceeds the any-T-age threshold) and Chr 9 at 107 Mb (*B*, 4.83 at T680); on high fat, proximal Chr 1 at 23.5 Mb (*D*, 4.82 at T20) and 37.6 Mb (*D*, 5.41 at T500). Combining each marker's T-ages (Cauchy combination; Methods), chow Chr 1 near 5 Mb passes over T20–T770 (−log10 *P* 5.30 against 4.68). In **c**–**f**, thresholds come from 200 permutations of strains within epochs and of F1 types within class (pending). Without the F1 mice (Extended Data), only chow Chr 1 near 5 Mb passes, at T770 (4.74 against 4.62). The shading is not corrected for the number of markers or T-ages; neighboring survivorships share most of their mice. The F1 mice were born in 2020–2022 and housed at TSRB; the F1 class covariates absorb any difference they share. Exploratory analyzes.

*Source: `qtl2_strainF1.R`, `qtl2_strainF1_perm.R` and `significant_loci_F1.R` (**a**,**b**); `qtl2_individualsF1.R` (**c**–**f**; thresholds from `QTL/qtl2_individualsF1_perm.R` in bxd_mortality_rates_2026); `cauchy_strainmeans.R`; panels by `figure_maps.R`.*

## Fig. 2 (earlier version, without the F1s; for Extended Data)

**Fig. 2 | Actuarial genetic maps of lifespan in female BXD mice on chow and high-fat diets.** Each map is a stack of genome scans, one per survivorship. A survivorship is the set of mice alive at a minimum age (T-age, left *y* axis). The maps run from T50 to T770 in 30-day steps (computed from T20), labeled every 60 days. Columns are the 21,056 markers of the BXD.8 genotypes on chromosomes 1–19 and X, placed by physical position.

Each scan is a linear mixed model (R/qtl2) with a polygenic effect whose kinship matrix leaves out the chromosome scanned (LOCO). In **a**,**b** the trait is the mean age at death of each strain with at least 3 mice alive at the T-age, with every strain weighted equally and no covariates. In **c**–**f** it is the age at death of individual mice, with covariates for vivarium and set-up age (180-day classes).

Shading gives −log10 *P* of the genotype term (likelihood-ratio test, 1 d.f.), in bins of 1–2, 2–3, 3–4, 4–5 and >5; white is below 1. Warm colors (yellow to red) mark markers at which mice with the *D* allele live longer, and cool colors (pale green to dark blue) markers at which mice with the *B* allele do. Ticks under each chromosome mark every 25 Mb, with longer ticks every 50 Mb. Numbers of strains (**a**,**b**) or mice (**c**–**f**) in each survivorship are given on the right *y* axis.

**a**, Strain means on chow diet (17% of calories from fat): 93 strains at T50 and 46 at T770.
**b**, Strain means on high-fat diet (60% of calories from fat): 87 strains at T50 and 14 at T770.
**c**, Individual mice on chow: *n* = 928 mice of 114 strains at T50 and 292 mice of 79 strains at T770.
**d**, Individual mice on high fat: *n* = 817 mice of 99 strains at T50 and 110 mice of 48 strains at T770.
**e**, Both diets, main effect of genotype, with diet as a further covariate: *n* = 1,745 mice of 119 strains at T50 and 402 mice of 86 strains at T770.
**f**, Genotype × diet interaction in the mice of **e**. Warm colors where the *D* − *B* difference in lifespan is larger on the high-fat diet, cool colors where it is smaller; −log10 *P* of the interaction term (likelihood-ratio test, 1 d.f.).

**Genome-wide significance (a,b).** Thresholds come from 1,000 permutations of strains, in which each strain took the whole actuarial series of another, so the nesting of survivorships was kept.
- **Thresholds, as −log10 *P* for a single T-age (for any marker at any T-age):**

  | Level | Chow | High fat |
  |---|---|---|
  | Significant, genome-wide *P* = 0.05 | 4.4 (5.5) | 4.4 (5.5) |
  | *P* = 0.10 | 4.1 (5.2) | 4.1 (5.0) |
  | Suggestive, *P* = 0.63 | 3.0 (4.0) | 3.0 (3.9) |

  The suggestive level is a result expected by chance about once per genome scan (Lander & Kruglyak 1995).
- **Thresholds hold only up to a T-age set by the number of strains:** T770 on chow (46 strains) and T620 on high fat (60 strains). The table's values cover T20 to those T-ages.
  - **In that range they are steady across T-ages:** 4.2–4.5 at 5% and 2.9–3.1 at the suggestive level.
  - **Beyond T620 on high fat (52 to 14 strains) they become unstable.** The 5% level runs from 3.5 to 5.9, and at T770 (14 strains) 70 of the 1,000 permutations tie at one maximum.
  - **So the high-fat rows from T650 to T770 are shown but not tested.**
- **Loci that pass the single-T-age threshold:**
  - on chow, Chr 1 at 5.0 Mb (*D* allele, 6.3 at T770, which also passes for any T-age) and Chr 9 at 107 Mb (*B*, 4.5 at T680);
  - on high fat, Chr 1 at 23.5 Mb (*D*, 5.0 at T20, passing from T20 to T110 and at T500) and Chr 2 at 157.5 Mb (*B*, 4.9 at T260, passing from T20 to T260).
- **Loci at the suggestive level only:**
  - on chow, Chr 2 at 157 Mb (*B*, 3.2 at T110), Chr 7 at 88 Mb (*B*, 3.2 at T740) and Chr 14 at 49 Mb (*B*, 3.4 at T380);
  - on high fat, Chr 3 at 64 Mb (*D*, 3.1 at T380), Chr 10 at 40 Mb (*D*, 3.8 at T500) and Chr 15 at 14 Mb (*B*, 3.5 at T50).

In **c**–**f**, kinship explains 33% of the variance in age at death at T20 on both diets (35% on chow, 38% on high fat), and 8–31% from T650 to T770. The strain means of **a**,**b** were also mapped weighted by the inverse of each mean's variance and by numbers of mice; both gave similar maps (*r* = 0.97–0.999 and 0.87–0.97 over markers at each T-age) with the same strongest loci (Methods). The shading shows −log10 *P* of each scan, not corrected for the number of markers or T-ages; neighboring survivorships share most of their mice. Exploratory analyzes.

*Source: `qtl2_strainmeans.R` (**a**,**b**; thresholds from `qtl2_permutations.R`, `output/qtl2_permutation_thresholds.txt`, and the loci that pass from `significant_loci.R`, `output/significant_loci.txt`; outputs `qtl2_*_means_*_Progressive.txt`, `qtl2_means_N.txt`, `qtl2_means_variance.txt`; the comparison of weightings in `compare_weights.R`, `output/compare_weights.txt` and `output/maps/Strain_means_three_weightings.svg`) and `qtl2_individuals.R` (**c**–**f**); panels `output/maps/Figure2_lifespan_six.svg`. Rob, 21 November 2026: keep equal weights in the figure, and note the other two weightings in the legend and Methods.*

**If Fig. 2 keeps Danny Arends's maps instead** (`output/maps/Figure2_lifespan.svg`, four panels), they read as follows:
- **a**,**b**, Strain means of age at death, for strains with at least 3 mice alive at the T-age, regressed on genotype (F test). Strains are weighted by their numbers of mice, with each strain's mean set-up age (180-day classes) as a covariate and no kinship. The right *y* axis gives the numbers of strains: 93 at T50 to 46 at T770 on chow, and 87 to 14 on high fat.
- **c**,**d**, Individual mice in a linear mixed model (lme4) with diet, vivarium and set-up age class as covariates and a random intercept for strain, but no kinship matrix. Likelihood-ratio tests; *n* = 1,745 mice at T50 to 402 at T770.

## Fig. 3: body weight and CTL

**Fig. 3 | Actuarial genetic maps of body weight, and of loci that modulate the correlation of body weight with lifespan, in female BXD mice.**

**a**,**b**, Body weight on chow (**a**) and high-fat (**b**) diets, from 80 to 740 days in 30-day steps. These are the ages at which more than 400 mice were weighed; the rows at 50 and 770 days are empty.
- Each mouse's weight is interpolated linearly between its weighings.
- Strain means, for strains with at least 3 mice weighed at the age, are regressed on genotype, with strains weighted by their numbers of mice (F test; no kinship).
- Warm colors: mice with the *D* allele are heavier; cool colors: mice with the *B* allele are heavier.
- The right *y* axis gives the numbers of strains: 43 at 80 days and 44 at 740 days on chow; 44 and 15 on high fat.

**c**,**d**, Correlated trait locus (CTL) maps on chow (**c**) and high-fat (**d**) diets.
- At each marker and age, body weight at that age and lifespan, each adjusted for vivarium and set-up age class, are rank-correlated separately in the mice that carry the *B* allele and in those that carry the *D* allele.
- Each mouse is weighted by the inverse of the number of mice of its strain.
- Shading gives −log10 *P* of the difference between the two correlations (χ² test on Fisher's *z*, 1 d.f.). Warm colors mark markers where the correlation is more positive in *D* carriers (*r*(*D*) > *r*(*B*)), and cool colors where it is more positive in *B* carriers.
- Within each genotype the correlation is mostly negative in younger survivorships and positive in older ones. Negative means heavier mice die sooner. It is negative up to about 530 days on chow and 290 days on high fat. So a warm color means a weaker negative correlation, or a stronger positive one, in *D* carriers.
- Mice entered a map if they were on their diet by the age and alive at it. Markers with fewer than 20 mice of either genotype are blank.

Rows, bins, colors, ticks and white shading are as in Fig. 2. *P* values are not corrected for the number of markers or ages. Exploratory analyzes.

*Source: `bodyweight.R`, `adjust.R` and `CTL.R` (run by `run_CTL.R`); Danny Arends's `QTL_BW_5Nov.xlsx`; the signs of the correlations from `ctl_signs.R` (`output/ctl_signs.txt`); panels `output/maps/Figure3_bodyweight_CTL.svg`.*

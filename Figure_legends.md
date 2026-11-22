# Figure legends: draft for the Nature Aging paper

*Draft by Claude (Opus 5.5) for Rob Williams, 21 November 2026. Not yet reviewed. The style follows Arends et al. 2026 (*Nature* 654:437): British (Oxford) spelling, "Fig. N |", and bold panel letters. Each number comes from the files named under its figure.*

*Before submission:*
- *A Nature legend title states the finding. These titles describe instead, because no genome-wide thresholds have been computed yet; permutations over strains would supply them.*
- *The small model notes under each panel belong in these legends and in Methods; drop them from the final panels.*

## Fig. 2: lifespan

**Fig. 2 | Actuarial genetic maps of lifespan in female BXD mice on chow and high-fat diets.** Each map is a stack of genome scans, one per survivorship. A survivorship is the set of mice alive at a minimum age (T-age, left *y* axis). The maps run from T50 to T770 in 30-day steps (computed from T20), labelled every 60 days. Columns are the 21,056 markers of the BXD.8 genotypes on chromosomes 1–19 and X, placed by physical position.

Each scan is a linear mixed model (R/qtl2) with a polygenic effect whose kinship matrix leaves out the chromosome scanned (LOCO). In **a**,**b** the trait is the mean age at death of each strain with at least 3 mice alive at the T-age, with every strain weighted equally and no covariates. In **c**–**f** it is the age at death of individual mice, with covariates for vivarium and set-up age (180-day classes).

Shading gives −log10 *P* of the genotype term (likelihood-ratio test, 1 d.f.), in bins of 1–2, 2–3, 3–4, 4–5 and >5; white is below 1. Warm colours (yellow to red) mark markers at which mice with the *D* allele live longer, and cool colours (pale green to dark blue) markers at which mice with the *B* allele do. Ticks under each chromosome mark every 25 Mb, with longer ticks every 50 Mb. Numbers of strains (**a**,**b**) or mice (**c**–**f**) in each survivorship are given on the right *y* axis.

**a**, Strain means on chow diet (17% of calories from fat): 93 strains at T50 and 46 at T770.
**b**, Strain means on high-fat diet (60% of calories from fat): 87 strains at T50 and 14 at T770.
**c**, Individual mice on chow: *n* = 928 mice of 114 strains at T50 and 292 mice of 79 strains at T770.
**d**, Individual mice on high fat: *n* = 817 mice of 99 strains at T50 and 110 mice of 48 strains at T770.
**e**, Both diets, main effect of genotype, with diet as a further covariate: *n* = 1,745 mice of 119 strains at T50 and 402 mice of 86 strains at T770.
**f**, Genotype × diet interaction in the mice of **e**. Warm colours where the *D* − *B* difference in lifespan is larger on the high-fat diet, cool colours where it is smaller; −log10 *P* of the interaction term (likelihood-ratio test, 1 d.f.).

In **c**–**f**, kinship explains 33% of the variance in age at death at T20 on both diets (35% on chow, 38% on high fat), and 8–31% from T650 to T770. Weighting the strain means of **a**,**b** by the inverse of their variance (between-strain variance not explained by kinship plus within-strain variance over *n*) gives nearly the same maps (*r* = 0.97–0.999 over markers at every T-age); weighting by numbers of mice moves some weak peaks (Extended Data Fig. N). *P* values are not corrected for the number of markers or T-ages, and neighbouring survivorships share most of their mice. Exploratory analyses.

*Source: `qtl2_strainmeans.R` (**a**,**b**; outputs `qtl2_*_means_*_Progressive.txt`, `qtl2_means_N.txt`, `qtl2_means_variance.txt`; the comparison of weightings in `compare_weights.R` and `output/maps/Strain_means_three_weightings.svg`) and `qtl2_individuals.R` (**c**–**f**); panels `output/maps/Figure2_lifespan_six.svg`.*

**If Fig. 2 keeps Danny Arends's maps instead** (`output/maps/Figure2_lifespan.svg`, four panels), they read as follows:
- **a**,**b**, Strain means of age at death, for strains with at least 3 mice alive at the T-age, regressed on genotype (F test). Strains are weighted by their numbers of mice, with each strain's mean set-up age (180-day classes) as a covariate and no kinship. The right *y* axis gives the numbers of strains: 93 at T50 to 46 at T770 on chow, and 87 to 14 on high fat.
- **c**,**d**, Individual mice in a linear mixed model (lme4) with diet, vivarium and set-up age class as covariates and a random intercept for strain, but no kinship matrix. Likelihood-ratio tests; *n* = 1,745 mice at T50 to 402 at T770.

## Fig. 3: body weight and CTL

**Fig. 3 | Actuarial genetic maps of body weight, and of loci that modulate the correlation of body weight with lifespan, in female BXD mice.**

**a**,**b**, Body weight on chow (**a**) and high-fat (**b**) diets, from 80 to 740 days in 30-day steps. These are the ages at which more than 400 mice were weighed; the rows at 50 and 770 days are empty.
- Each mouse's weight is interpolated linearly between its weighings.
- Strain means, for strains with at least 3 mice weighed at the age, are regressed on genotype, with strains weighted by their numbers of mice (F test; no kinship).
- Warm colours: mice with the *D* allele are heavier; cool colours: mice with the *B* allele are heavier.
- The right *y* axis gives the numbers of strains: 43 at 80 days and 44 at 740 days on chow; 44 and 15 on high fat.

**c**,**d**, Correlated trait locus (CTL) maps on chow (**c**) and high-fat (**d**) diets.
- At each marker and age, body weight at that age and lifespan, each adjusted for vivarium and set-up age class, are rank-correlated separately in the mice that carry the *B* allele and in those that carry the *D* allele.
- Each mouse is weighted by the inverse of the number of mice of its strain.
- Shading gives −log10 *P* of the difference between the two correlations (χ² test on Fisher's *z*, 1 d.f.). Warm colours mark markers where the correlation is more positive in *D* carriers (*r*(*D*) > *r*(*B*)), and cool colours where it is more positive in *B* carriers.
- Within each genotype the correlation is mostly negative in younger survivorships and positive in older ones. Negative means heavier mice die sooner. It is negative up to about 530 days on chow and 290 days on high fat. So a warm colour means a weaker negative correlation, or a stronger positive one, in *D* carriers.
- Mice entered a map if they were on their diet by the age and alive at it. Markers with fewer than 20 mice of either genotype are blank.

Rows, bins, colours, ticks and white shading are as in Fig. 2. *P* values are not corrected for the number of markers or ages. Exploratory analyses.

*Source: `bodyweight.R`, `adjust.R` and `CTL.R` (run by `run_CTL.R`); Danny Arends's `QTL_BW_5Nov.xlsx`; the signs of the correlations from `ctl_signs.R` (`output/ctl_signs.txt`); panels `output/maps/Figure3_bodyweight_CTL.svg`.*

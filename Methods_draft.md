# Methods: draft paragraphs for the Nature Aging paper

*Draft by Claude (Opus 5.5) for Rob Williams, 21 November 2026. Not yet reviewed. The numbers come from the outputs named at the end.*

## Actuarial mapping of lifespan

**Survivorships.** For each diet and each T-age from 20 to 770 days in 30-day steps, the survivorship is the set of female BXD mice alive at the T-age, that is, with an age at death at or above it. Maps are shown from T50.

**Genotypes.** We used the 21,056 markers of the BXD.8 genotypes (GeneNetwork; 603 on the X), with heterozygous and unknown calls treated as missing. The BXD strains were treated as a sib-mated recombinant inbred panel in R/qtl2 (version 0.46; Broman et al. 2019), with genotype probabilities from `calc_genoprob` (error probability 0.002). The X chromosome was treated like an autosome.

**Kinship.** Kinship was the proportion of alleles shared, computed from the strains' genotype probabilities, leaving out the chromosome scanned (LOCO).

**Strain means.** The trait was the mean age at death of the mice of each strain in the survivorship, for strains with at least 3 such mice.
- Every scan was a linear mixed model with a polygenic effect and LOCO kinship (`scan1`, REML), with no covariates.
- Strain means were weighted equally.
- We also mapped them with two other weightings:
  - by numbers of mice;
  - by the inverse of the variance of each strain mean, 1/(σ²ₛ + σ²_w/*n*). Here σ²_w is the variance among mice of a strain, pooled over strains, and σ²ₛ is the variance between strains not explained by kinship. σ²ₛ was estimated by REML for the strain means, with covariance σ²ₖ**K** + diag(σ²ₛ + σ²_w/*n*).
- The three weightings gave similar maps. Over markers, at each T-age, equal weights correlated with inverse-variance weights at *r* = 0.97–0.999 and with weights by numbers of mice at *r* = 0.87–0.97. The strongest loci were the same under all three.

**Individual mice.** The trait was the age at death of each mouse, and each mouse took its strain's genotype probabilities and kinship.
- Covariates were vivarium and set-up age in 180-day classes, with diet as a further covariate when both diets were mapped together.
- The genotype × diet interaction was tested by the difference between the log-likelihoods of models with and without the interaction term.

**Scores.** Each scan's LOD was converted to *P* with the χ² distribution on 1 d.f. Maps show −log10 *P*, signed by the direction of the *D* − *B* effect from `scan1coef` with the same model.

*Sources: `qtl2_cross.R`, `qtl2_strainmeans.R`, `qtl2_individuals.R` and `compare_weights.R` (outputs `compare_weights.txt` and `qtl2_means_variance.txt`) in robwwilliams/BxDLifeSpan.*

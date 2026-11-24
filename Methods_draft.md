# Methods: draft paragraphs for the Nature Aging paper

*Draft by Claude (Opus 5.5) for Rob Williams, 21 November 2026; F1 hybrids, the Cauchy combination and two-locus scans added 23 November 2026. Not yet reviewed. The numbers come from the outputs named at the end.*

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

**Genome-wide thresholds.** For the strain-mean maps, thresholds came from 1,000 permutations of strains, drawn from a fixed seed.
- **In each permutation,** every strain took the strain means of another strain at all T-ages, so the nesting of the survivorships was kept and only the link between genotype and lifespan was broken. Strains kept their genotypes and kinship.
- **The largest −log10 *P*** over markers was recorded at each T-age, and over markers and T-ages together.
- **Stable range.** Thresholds were used only from T20 to the last T-age before the 5% single-T-age threshold first moved more than 0.3 from its median over T20–T500. That is T770 on chow (46 strains) and T620 on high fat (60 strains). Beyond T620 on high fat, with 52 to 14 strains, the thresholds were unstable, and at T770 (14 strains) 70 of 1,000 permutations tied at one maximum.
- **The 95th percentiles** gave a 5% threshold for a single T-age (median 4.4 on both diets, 4.2–4.5 across the stable range), and for any marker at any T-age in that range (5.5 on both diets).
- **The 90th and 37th percentiles** gave the 10% level (4.1 for one T-age; 5.2 on chow and 5.0 on high fat for any T-age) and the suggestive level of Lander & Kruglyak (1995), genome-wide *P* = 0.63 (3.0 for one T-age; 4.0 and 3.9 for any T-age).

**BXD × BXD F1 hybrids.** Thirteen F1 types (6 short × short and 7 long × long crosses, chosen for recombinations over proximal Chr 1) were mapped with the strains; 308 female F1 mice died naturally on chow or high fat, all born in 2020–2022 and housed at TSRB.
- **Genotypes.** Each F1 type took its genotype probabilities from both parents: P(BB) = p_B1 p_B2, P(BD) = p_B1 p_D2 + p_D1 p_B2, P(DD) = p_D1 p_D2. Additive scans used the *D* allele dosage. Kinship was allele sharing over strains and F1 types together, leaving out the chromosome scanned.
- **Strain and F1 means.** Each F1 type was one unit, like a strain: the mean age at death of its mice alive at the T-age, with at least 3 such mice, every unit weighted equally. F1 class (inbred, short × short, long × long) was a covariate, so that the average hybrid vigour of a class was not read as a genotype effect. The dominance deviation, BD − (BB + DD)/2, was tested by the difference between a three-genotype scan and the additive scan (1 d.f.).
- **Individual mice.** BXD and F1 mice were mapped together with the covariates of the BXD maps, plus F1 class and, for the maps of both diets, F1 class × diet. The diet term for F1 class keeps the difference in hybrid vigour between diets (about 196 days on chow and 146 on high fat for short × short F1s) out of the genotype × diet test. The F1 workbook records no vivarium; with the F1 class covariates in the model, any coding of it gives the same fit.

**Thresholds for the strain and F1 maps.** 1,000 permutations: each strain took the whole actuarial series of another strain, and each F1 type that of another F1 type of its class. The stable range of T-ages, by the rule above, was T20–T740 on chow and T20–T590 on high fat for the additive maps (T20–T500 for dominance). The 5% thresholds were 4.66 on chow and 4.59 on high fat for one T-age (5.83 and 5.46 for any T-age), and 1.15 and 1.91 for dominance; only the 13 F1 types inform dominance, so its −log10 *P* cannot be read on the scale of a test on 1 d.f. These permutations break each F1's resemblance to its parents, so for dominance they may set thresholds too low.

**Thresholds for the individual-mouse maps** came from 200 permutations in which each strain took the genotypes and kinship of another strain of its stratum (Taylor's strains, BXD1–BXD42, or the newer ones) and, with F1s, each F1 type those of another F1 type of its class. Every mouse kept its lifespan and covariates. Each permutation took about 25 minutes per core; they were run on a lab workstation in chunks with fixed seeds.

**One P per marker over T-ages.** The P values of a marker's T-ages were combined with the Cauchy combination test (Liu & Xie 2020), as in Arends et al. 2026, over the stable range of T-ages, with exact P values on 1 d.f. (none set to 1). Genome-wide thresholds for the largest combined −log10 *P* came from 1,000 permutations, drawn as above. At 5% they were 4.44 (chow) and 4.41 (high fat) for strain means, and 4.82 and 4.73 for strain and F1 means.

**Support intervals and two loci.** Support intervals were the runs of neighbouring markers whose LOD came within 1.5 of the peak's. Whether proximal Chr 1 on high fat holds one locus or two was tested at each T-age from T20 to T620 by fitting every pair of positions in 0–50 Mb in the same mixed model (heritability estimated without loci and held fixed; generalized least squares) and comparing the best pair with the best single locus; P came from 1,000 parametric bootstraps under the best single locus, each put through the whole search.

*Sources: `qtl2_cross.R`, `qtl2_strainmeans.R`, `qtl2_individuals.R`, `compare_weights.R` (outputs `compare_weights.txt` and `qtl2_means_variance.txt`), `qtl2_permutations.R` (`qtl2_permutation_thresholds.txt`) and `significant_loci.R` (`significant_loci.txt`) in robwwilliams/BxDLifeSpan. With F1s: `strainF1_setup.R`, `qtl2_strainF1.R`, `qtl2_strainF1_perm.R` and `significant_loci_F1.R`, `qtl2_individualsF1.R`, `cauchy_strainmeans.R`, `chr1_zoom.R`, `chr1_f1_dropone.R` and `chr1_twolocus.R`; the permutations of the individual-mouse maps are `QTL/qtl2_individuals_perm.R` and `QTL/qtl2_individualsF1_perm.R` in robwwilliams/bxd_mortality_rates_2026.*

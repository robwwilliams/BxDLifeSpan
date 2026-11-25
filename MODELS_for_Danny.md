# The models in this fork, for Danny

*Written by Claude (Opus 5.5) for Rob Williams and Danny Arends, 24 November 2026. Rob: "As long as any new code gives results that are generally harmonious with less complex models, then all is good." Danny's own files are unchanged in this fork.*

## 1. What is new, in the language of `lm`

| Map | Script | Roughly, in lme4 terms |
|---|---|---|
| Strain means | `qtl2_strainmeans.R` | `lm(mean ~ D)` + a polygenic term |
| Strain and F1 means | `qtl2_strainF1.R` | `lm(mean ~ F1class + D)` + a polygenic term |
| Individual mice | `qtl2_individuals.R` | `lmer(age ~ vivarium + setup + D + (1 \| strain))`, with the kinship matrix in place of the strain intercept |
| Individual mice and F1s | `qtl2_individualsF1.R` | the same, plus `F1class` (and `F1class:diet` for both diets) |
| Genotype × diet | same scripts | the likelihood ratio of models with and without `D:diet` |

- **`D`** is the share of a unit's alleles that are *D*: 0 (*BB*), 0.5 (*BD*, only an F1) or 1 (*DD*). Its coefficient is *DD* minus *BB*, in days.
- **"Polygenic term":** a random effect whose covariance is the kinship matrix, estimated by REML, leaving out the chromosome scanned (LOCO). It replaces the random strain intercept: strains that share much of their genome get correlated intercepts.
- **R/qtl2 doubles the kinship internally** (variance = *h* · 2*K* + (1 − *h*) *I*). This matters only if you recompute a LOD by hand.

## 2. The F1 part, which is the unfamiliar bit

**Read [`strainF1_setup.R`](strainF1_setup.R) first;** it is 70 lines, and everything else uses it.
1. **An F1 type has no genotypes of its own.** It takes its two parents' probabilities, marker by marker: P(*BD*) = P(parent 1 *B*) P(parent 2 *D*) + P(parent 1 *D*) P(parent 2 *B*). So it is *BD* exactly where its parents differ.
2. **Each F1 type is one row** of the strain-mean maps, with the same weight as a strain.
3. **F1 class** (inbred, short × short, long × long) is a covariate. Without it, the F1s' hybrid vigor would be read as an effect of every marker at which F1s differ from strains.
4. **Kinship is computed over strains and F1s together,** so an F1 is correlated with its parents and with other F1s that share a parent.
5. **Dominance** is the LOD of a three-genotype scan minus the additive LOD. Only the 13 F1 types inform it.

**Checks on this part:** [`QTL/test_alignment.R`](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/QTL/test_alignment.R) in the other repository tests that mice, units, genotypes and kinship line up, with deliberate off-by-one errors that must be caught. It covers the F1 permutation script only.

## 3. Do they agree with simpler models? ([`harmony_F1_models.R`](harmony_F1_models.R))

Correlation over markers of the signed −log10 P, at T20, T380, T500 and T620:

| Comparison | Chow | High fat |
|---|---|---|
| Strains alone: mixed model against plain `lm` | 0.98–1.00 | 0.95–1.00 |
| Strains + F1s: mixed model against plain `lm` with F1 class | 0.95–0.98 | **0.84–0.96** |
| Mixed model: with F1s against without | 0.96–0.99 | 0.95–0.98 |
| Individual mice with F1s: R/qtl2 against `lmer` with a strain intercept (300 markers; T20 and T500) | 0.81–0.87 | 0.86 |

- **Harmonious:** three to five of the five strongest chromosomes are the same in each pair.
- **The largest gap is on high fat with F1s at T20 and T380 (0.88 and 0.84).** There the kinship term matters most: the F1s share parents, and plain `lm` treats them as independent.

## 4. Your CTL map, run both ways ([`ctl_permutations.R`](ctl_permutations.R))

**Your way is reproduced exactly:** 1,052,572 values at 25 T-ages and two diets differ from your `CTL_<diet>_LODmatrix.txt` only by rounding (at most 0.005), and the same cells are missing.

**Then the corrections of Arends et al. 2026, and permutations of strains** (1,000; each strain takes another strain's genotypes; every mouse keeps its weight, lifespan and strain):

| | Chow | High fat |
|---|---|---|
| Inflation (median χ² / expected), real data | 2.75 | 2.51 |
| Inflation, strains permuted (no CTL at all) | 2.55 (90% range 2.20–2.98) | 2.24 (1.95–2.55) |
| Inflation, (weight, lifespan) pairs shuffled among mice | 1.77 | 1.61 |
| B&H, 5% FDR over each T-age's markers: T-ages with calls | 21 of 25 | 23 of 25 |
| Permuted data with at least one B&H call | 100% | 100% |
| Permuted data with at least one Bonferroni call (all markers and T-ages) | 100% | 99% |
| Largest −log10 P, any T-age | 10.2 | 9.4 |
| 5% permutation threshold for it | 16.8 | 12.3 |
| Genome-wide P of the largest | 0.88 | 0.50 |

- **The BXD map looks as strong after the genotypes are scrambled among strains as before.** B&H and Bonferroni cannot see this, because both start from P values that count mice.
- **Two causes:** the weights 1 / (mice per strain) leave fewer effective mice than the n in Fisher's z (inflation 1.8 and 1.6 with no strain structure at all), and the mice of a strain share weight, lifespan and genotype (the rest).
- **One T-age is below 5% on its own:** high fat T470, Chr 18 at 18.7 Mb, 9.4 against 8.8 (P = 0.030). With 50 T-age × diet scans it is a lead at best, not a finding.
- **This does not touch Arends et al. 2026.** In UM-HET3 every mouse is genetically unique, so n is right there. The rule does not carry over to strains with many identical mice.
- **The test can detect a real CTL, but only a large one.** Planted at Chr 6, 140 Mb (T200 to T500), a difference of 50 or 100 days of lifespan per SD of weight in *D* mice was missed, and 200 days per SD was found only in part: on high fat at T200 and T290 (P = 0.005 and 0.010) but not at T380 and T500, and on chow at T500 (0.035; 0.11–0.53 at the other three T-ages). Our 119 strains carry little power for CTLs. (One marker and 200 permutations per size, so this is a rough power check.)
- **The permutation swaps strains within two strata (Taylor's and newer) and takes no account of kinship** beyond that.
- **Within strains only** (`ctl_permutations.R 1000 12 within`; a different question, not your statistic): chow real 6.3 against a threshold of 8.7 (P = 0.68); high fat 8.7 against 11.9 (P = 0.37).

**Suggestions for the paper:** show your CTL maps in Extended Data as exploratory, say that they were tested by permuting strains, and leave the Results claim out; or look for a design with more power, such as a model of weight × genotype with a strain intercept and slope (not yet tried).

**Not done (stretch goals, Rob, 24 November 2026):** CTL × diet, and epistasis (Vita1a × Vita6b, 9a and 11a in BXD females). For CTL × diet, power will be below that of either diet alone.

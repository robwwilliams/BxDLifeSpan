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

### The code cleaned up, at your request (26 November 2026)

`ctl_permutations.R` is replaced by four short files. Each has one path and takes no command-line input.

| File | What it does |
|---|---|
| [`ctl_inputs.R`](ctl_inputs.R) | Builds the two input tables (`output/ctl_scan/mice.txt`, `genotypes.txt`) through your `shared.R` and `adjust.R`, unchanged |
| [`ctl_functions.R`](ctl_functions.R) | One function per step, numbered as below |
| [`ctl_scan.R`](ctl_scan.R) | The settings, then the steps in order; writes `output/ctl_scan/` |
| [`ctl_tests.R`](ctl_tests.R) | Checks every function against your code or a committed result; writes `output/ctl_scan/tests.txt` |

**The approach, in six steps:**
1. **Inputs.** Lifespan and weight of the BXD mice on chow or high fat, adjusted by your `adjust.R` within each diet: `~ vivarium + set-up age class` (0–180, 180–360, … days), residual + mean, rounded to a day. BXD mice only; no F1s.
2. **Scans.** A scan is one diet, one survivorship age T and one weighing age W. Its mice were on the diet by W and alive at T.
   - Your design (`actuarial`) weighs at T.
   - The grid designs (`W80` to `W740`) weigh at one age and follow T from W to 770 days.
3. **Map.** Your test at each marker and scan: the weighted Spearman correlations in the B and the D mice, compared by Fisher's z and your chi-square on 1 d.f. It gives a signed −log10 P.
4. **Two statistics per design and diet.** They are separate tests and are not expected to be equal:
   - **max:** the largest −log10 P over all markers and T-ages;
   - **cauchy:** at each marker, its T-ages' P values joined by your `cauchycombination()`, then the largest over markers.
5. **Permutations.** Each strain takes the genotypes of another strain of its stratum (BXD1–42, or the newer strains), and the whole map is rescanned.
   - One assignment serves every scan, diet and design.
   - Mice keep their weights, lifespans and strain mates.
6. **Verdict.** Each statistic of the real map is judged against the same statistic of the 1,000 permuted maps: a 5% threshold and a genome-wide P.

**Your six points:**
1. **F1s are dropped.**
2. **The linear adjustments** are your `adjust.R`, called in `ctl_inputs.R`, with the models written out in its header.
3. **Structure.** The steps are named functions, and `ctl_scan.R` is the settings plus about 40 lines that call them.
4. **The `atanh()`.** It is gone: `ctlTest()` uses your own expression, `0.5 * log((1 + r) / (1 - r))`. The tests:
   - compare it with your loop body from `CTL.R` on 5,000 random pairs (identical);
   - r = 0 gives z = 0;
   - |r| = 1 gives NA (`CTL.R` gives NaN);
   - P below about 1e-308 gives infinity in both versions.
5. **The approach** is the six steps above.
6. **No command-line inputs.** The settings are at the top of `ctl_scan.R`. `run_info.txt` records them, with the md5 of the input tables and of the code.

**The tests** (`output/ctl_scan/tests.txt`) all pass:
- **The C++ correlation** equals your `cov.wt(apply(mat, 2, rank), wt = weights, cor = TRUE)` line, both on made-up data with ties and at 100 real markers (difference below 1e-15).
- **The real map** equals your `CTL_CD_LODmatrix.txt` and `CTL_HF_LODmatrix.txt` to the 2 decimals you saved, with the same cells missing.
- **Six maps of the committed run (the real one and five permutations)** come out the same:
  - the largest −log10 P, the Bonferroni and B&H calls and the inflation agree to 4 decimals;
  - the Cauchy maximum agrees within 0.013.

**One fault found in the old code.** `ctl_permutations.R` computed the Cauchy quantile as `tan((0.5 - p) * pi)`. This loses digits for P between about 1e-15 and 1e-12: at P = 1e-15 the quantile is off by 0.9%. Your `qcauchy()`, now used, is exact.
- The real map is unaffected.
- One permuted map of the six checked moves by 0.012.
- The 1,000-permutation rerun will give the exact thresholds.

**Why 1,000 permutations took minutes and `CTL.R` takes hours.** Every permutation does rescan every marker, T-age and diet. Three things make it fast:
- the correlations run in C++;
- markers with the same strain pattern are computed once;
- the permutations run in parallel.

`CTL.R` also writes and re-reads a table at every T-age.

**How the old code handled the F1 mice, and why they are now dropped.**
- **Genotype.** An F1 type counted as B (or D) at a marker where both its parents are B (or D), because there it is homozygous.
- **Heterozygous markers.** Where the parents differ, the F1 is B/D. It was left out of both groups at that marker, as `CTL.R` leaves out an H or unknown call in a BXD strain.
- **No H group.** There was no third group and no model of heterozygotes. About half the markers are homozygous for each F1 type (0.45 to 0.58; mean 0.50, weighted by mice), so each F1 type entered the B or D group at about half the markers and was absent at the rest.
- **Centering.** Lifespan and weight of the F1 mice were centered on their class (short × short, long × long) within diet, so that hybrid vigor would not be read as a CTL.
- **Permutations.** F1 types swapped genotypes only within their class.
- **Why they are dropped.** The F1 mice are heterozygous across about half the genome, at markers other than the one tested, and centering by class is a rough correction for that. This is the reason the clean code leaves them out.

**Still to run.** The 1,000 permutations of every design, BXD only, take about an hour on 8 cores and about 4 GB. Until they are back, the grid results below include the F1s.

**Your way is reproduced exactly:** 1,052,572 values at 25 T-ages and two diets differ from your `CTL_<diet>_LODmatrix.txt` only by rounding (at most 0.005), and the same cells are missing.

**Then the corrections of Arends et al. 2026, and permutations of strains** (1,000; each strain takes another strain's genotypes; every mouse keeps its weight, lifespan and strain):

| | Chow | High fat |
|---|---|---|
| Inflation (median χ² / expected), real data | 2.75 | 2.51 |
| Inflation, strains permuted (no CTL at all) | 2.55 (90% range 2.20–2.98) | 2.24 (1.95–2.55) |
| Inflation, any mouse permuted with any mouse (naive) | 1.72 | 1.52 |
| B&H, 5% FDR over each T-age's markers: T-ages with calls | 21 of 25 | 23 of 25 |
| Permuted data with at least one B&H call | 100% | 100% |
| Permuted data with at least one Bonferroni call (all markers and T-ages) | 100% | 99% |
| Largest −log10 P, any T-age | 10.2 | 9.4 |
| 5% permutation threshold for it | 16.8 | 12.3 |
| Genome-wide P of the largest | 0.88 | 0.50 |

- **The BXD map looks as strong after the genotypes are scrambled among strains as before.** B&H and Bonferroni cannot see this, because both start from P values that count mice.
- **Two causes:** the weights 1 / (mice per strain) leave fewer effective mice than the n in Fisher's z (inflation 1.7 and 1.5 with no strain structure at all), and the mice of a strain share weight, lifespan and genotype (the rest).
- **A naive permutation, any mouse with any mouse, calls both maps significant** (P = 0.023 on chow, 0.019 on high fat), because it breaks up the strains. Its 5% threshold is exceeded by 96% (chow) and 79% (high fat) of maps with scrambled strains. Why, and how to avoid it: [Hierarchical_permutation_explained.md](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/notes/Hierarchical_permutation_explained.md).
- **With the F1 mice added (25 November 2026; `F1` as a further argument), the answer is the same:** largest 10.0 on chow (threshold 16.5, P = 0.92) and 10.3 on high fat (12.9, P = 0.35). An F1 type counts as *B* or *D* only where its parents agree, is centered on its class within diet, and is permuted within class. Its weights have not been through our full QC.
- **One T-age is below 5% on its own:** high fat T470, Chr 18 at 18.7 Mb, 9.4 against 8.8 (P = 0.030). With 50 T-age × diet scans it is a lead at best, not a finding. (With the F1s it is 10.3 against 8.9, P = 0.009.)
- **This does not touch Arends et al. 2026.** In UM-HET3 every mouse is genetically unique, so n is right there. The rule does not carry over to strains with many identical mice.
- **The test can detect a real CTL, but only a large one.** Planted at Chr 6, 140 Mb (T200 to T500), a difference of 50 or 100 days of lifespan per SD of weight in *D* mice was missed, and 200 days per SD was found only in part: on high fat at T200 and T290 (P = 0.005 and 0.010) but not at T380 and T500, and on chow at T500 (0.035; 0.11–0.53 at the other three T-ages). Our 119 strains carry little power for CTLs. (One marker and 200 permutations per size, so this is a rough power check.)
- **The permutation swaps strains within two strata (Taylor's and newer) and takes no account of kinship** beyond that.
- **Within strains only** (`ctl_permutations.R 1000 12 within`; a different question, not your statistic): chow real 6.3 against a threshold of 8.7 (P = 0.68); high fat 8.7 against 11.9 (P = 0.37).

### The design of Arends et al. 2026: weight at one age, survivorship rising (25 November 2026)

Rob: "did you take body weight at 3 months and then compute lifespan:weight correlations of the survivorship up toward 800 days?" The runs above did not: as in `CTL.R`, weight is taken at the T-age itself, so weight and survivorship move together. The other design holds weight at 90, 180, 365 or 545 days and steps only the survivorship, every 30 days (`weight_lifespan_fixed.R`, and `ctl_permutations.R 1000 4 fixed 90`).

- **The change of sign goes away.** With weight held at 90 or 180 days, the pooled correlation is negative at nearly every T-age, up to T770: heavier young mice die sooner, whichever survivors are asked. No T-age on either diet, at any of the four weight ages, has an interval wholly above zero. So the change of sign in section 4 comes from the age of the weight, late weights meaning something else (wasting), and not from the survivorship.
- **It is established on high fat.** Weight at 90 days, within strains: −0.15 at T90, −0.23 (95% CI −0.36 to −0.11) at T380, and no longer detectable from T500 on, where survivors are few. Weight at 180 days, pooled: −0.15 (−0.26 to −0.03) at T180 to −0.24 (−0.36 to −0.11) at T590.
- **On chow it is weaker:** mostly negative, and pooled intervals exclude zero at 7 of 21 T-ages for weight at 180 days and at one or none for the others.
- **Caution:** only 309 chow and 329 high-fat mice have a weight at 90 days, half the numbers at 180 days, and few survive to late T-ages (HF at T680: 84).
- **CTL maps, same design, 1,000 permutations of strains:** nothing passes genome-wide at any weight age (largest P 0.21 on high fat, weight at 180 days; 0.33 to 0.70 otherwise). Four of 124 single scans fall below 5% on their own, about 6 expected by chance. B&H again calls markers in 100% of the permuted maps.

### The full scan (25 November 2026)

Rob: "the full actuarial scan across all ages at which animals were weighed (overkill)". Weight at every 30 days from 80 to 740, against every later survivorship to T770; BXD and F1 mice; 1,000 hierarchical permutations, the same at every weight age (`ctl_permutations.R 1000 3 fixed <W> F1`, `ctl_grid_summary.R`, `ctl_grid_figure.R`; figure `output/maps/CTL_grid.svg`).

- **No CTL passes.** Over all 595 scans the largest −log10 P is 12.4 on chow (5% threshold 17.6, P = 0.68) and 11.6 on high fat (13.8, P = 0.31).
- **Single scans below 5%: 35 of 595, about 30 expected by chance** (11 on chow, 24 on high fat). Neighboring scans share mice and weights, so the colored cells come in runs; the count is not unusual on either diet (P = 0.49 and 0.19; in the permutations the median count is 10 and 12, and the 95th percentile 46 and 39).
- **Strongest leads:** high fat Chr 6 at 17.6 Mb, weight at 200–290 days, T620 (P = 0.007–0.019 in single scans), and chow Chr 11 at 111 Mb at T770, weight at 470–530 days (P = 0.004–0.014).

### Danny's questions: the actuarial scan, and why Bonferroni is so far off (25 November 2026)

Danny: "I wouldnt have exected the permutations to be so far off the BF adjustment... Did you permute the whole actuarial T(age) scan? If so, how did you deal with the correlation across ages ? did you use a Cauchi Correction ?"

- **Yes, the whole scan.** One permutation moves the genotypes of whole strains, once, and then reruns every T-age and both diets with that assignment. Mice keep their weights and lifespans, so the correlation across ages (and among linked markers) is the same in every permutation as in the data.
- **The correlation is handled by the maximum.** Each permutation keeps its largest −log10 P over all markers and T-ages. No independence is assumed and no number of effective tests is needed.
- **We had not used the Cauchy combination. It is now added, with the same permutations, and the answer is the same** (1,000 permutations; exploratory, after the results above). The Cauchy code gives the same numbers as `cauchycombination()` in your `Cauchy.R` (difference 4e-15), on the 1 d.f. P values of the CTL test.

| | Chow | High fat |
|---|---|---|
| Largest Cauchy −log10 P over markers | 8.8 (Chr 11, 99.7 Mb) | 8.0 (Chr 18, 18.7 Mb) |
| Bonferroni threshold for it, 5% over 21,056 markers | 5.6 | 5.6 |
| 5% permutation threshold, whole strains | 15.4 | 11.1 |
| Genome-wide P, whole strains (any T-age, as above) | 0.89 (0.88) | 0.53 (0.50) |
| Scrambled strains with at least one Bonferroni call on the Cauchy P | 100% | 100% |
| Same, B&H | 100% | 100% |

- **Cauchy cannot repair this.** It allows for the correlation among T-ages, but it takes each P value as correct, as Bonferroni does. Here each one is too small before any correction: the inflation is 2.55 and 2.24 with no CTL at all. The problem is the number of independent genomes, not the number of tests.
- **A control: where the unit is right, Bonferroni is right.** Shuffle any mouse with any mouse and give every mouse weight 1 (mode `micewt1`):

| No CTL at all; inflation (median χ² / expected) | Chow | High fat |
|---|---|---|
| Mice shuffled, weights 1 | 1.02 | 0.99 |
| Mice shuffled, weights 1 / (mice of the strain) | 1.72 | 1.52 |
| Whole strains shuffled, weights 1 / (mice of the strain): the real design | 2.55 | 2.24 |

  In the first row, Bonferroni over all markers and T-ages calls something in 9 (chow) and 3 (high fat) of 1,000 maps, Bonferroni on the Cauchy P in 10 and 6, and B&H on the Cauchy P in 52 and 38: none above its nominal 5% by more than chance, and Bonferroni is conservative, as expected for correlated tests. So Bonferroni is off for two separate reasons, which the next two rows add one at a time.
  1. **Weights.** At T200 on chow, 623 mice with these weights carry the information of 359 equally weighted ones (Kish), yet Fisher's z uses 623. The ratio, 1.74, is the chow inflation in row 2.
  2. **Shared genomes.** The 623 mice are 95 strains, and the intraclass correlation of mice within strains is 0.29 for lifespan and 0.52 for weight.

- **With the F1 mice:** Cauchy P = 0.92 on chow and 0.37 on high fat (the maximum: 0.92 and 0.35).
- **The full grid was not redone with Cauchy.** The rerun was stopped on 25 November 2026 after 17 of 23 weighing ages (80–560 days), at Rob's word that the CTL work is closed ("Nothing is significant and this work is only going into a paper as a note that sample size is inadequate"). The grid results above (595 scans, chow P = 0.68, high fat P = 0.31) are complete and are what that note rests on.

**Suggestions for the paper:** show your CTL maps in Extended Data as exploratory, say that they were tested by permuting strains, and leave the Results claim out; or look for a design with more power, such as a model of weight × genotype with a strain intercept and slope (not yet tried).

**Not done (stretch goals, Rob, 24 November 2026):** CTL × diet, and epistasis (Vita1a × Vita6b, 9a and 11a in BXD females). For CTL × diet, power will be below that of either diet alone.

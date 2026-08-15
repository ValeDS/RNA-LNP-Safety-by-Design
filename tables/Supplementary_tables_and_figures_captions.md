# Supplementary table and figure captions

## Supplementary tables

**Table S1. Complete synthetic response-generation model.** Equations, coefficients, additive Gaussian error distributions, nominal error variances and observed non-negativity truncation counts for the five synthetic endpoints. Adaptive activation is calculated using the realized delivery value, and cytokine burden is calculated using the realized innate-activation value; the model is therefore sequential. All quantities are synthetic and assumption-defined rather than experimentally calibrated.

**Table S2. Endpoint normalization and Safety-by-Design score definition.** Endpoint-specific minima, maxima and ranges estimated once from the primary 10,000-record dataset and used for min–max normalization. The composite score assigns weights of 0.30, 0.20, 0.20, 0.15 and 0.15 to delivery, adaptive activation, innate safety, cytokine safety and off-target safety, respectively.

**Table S3. Random Forest grouped cross-validation performance.** Pooled out-of-fold metrics from five-fold cross-validation grouped by formulation and, as a sensitivity analysis, by virtual patient. Formulation-grouped validation prevents the same formulation from appearing in training and validation and is the primary test. Performance quantifies recovery of the prespecified synthetic response surface, not experimental or clinical predictive validity.

**Table S4. Safety-by-Design ranking sensitivity to component weights.** Rank range and frequency of top-1, top-10 and top-20 membership for all 500 formulations across 3,876 weight combinations on a 0.05-step simplex, with every component constrained to a minimum weight of 0.05. The original leading candidate, LNP_0462, remained in the top 10 under 99.85% and in the top 20 under 100% of scenarios.

**Table S5. Balanced full-factorial inflammatory-risk analysis.** Part A reports N, mean, SD, median, first and third quartiles, minimum and maximum for all five endpoints and the Safety-by-Design score in each risk stratum across the complete 200-patient × 500-formulation design. Part B reports the complete rankings and Pareto membership of all 500 formulations in each stratum. Every formulation was evaluated in 109 Low-, 78 Intermediate- and 13 High-risk virtual patients.

**Table S6. Complete primary Pareto-optimal candidate set.** Seven formulation parameters, five endpoint means, mean and standard deviation of the Safety-by-Design score and record counts for all 21 independently verified non-dominated candidates. Rank denotes post hoc ordering by the composite score and is conceptually distinct from Pareto membership.

**Table S7. BPCI exact shuffled-label null and component-weight sensitivity.** Observed Biological Program Concordance Index values, exact null summaries and upper-tail p-values obtained by enumerating all 120 permutations of five program labels. Minimum, median and maximum BPCI values summarize 120 component-weight combinations with each component weight constrained to at least 0.10.

**Table S8. BPCI sample sizes and paired-donor bootstrap uncertainty.** Donor and cell counts by post-vaccination day and median BPCI with 95% percentile intervals from 2,000 donor-resampling bootstrap replicates. Each post-vaccination donor profile was compared with that donor's Day-0 baseline. The original estimator uses the overall Day-0 mean and is shown separately.

**Table S9. Synthetic multi-omics residual-noise audit after seed correction.** Empirical residual means, standard deviations and variances for the cytokine and immune-cell layers regenerated with independent reproducible seeds. Maximum residual correlations within and between layers confirm elimination of the perfect cross-layer correlations caused by reuse of the original random seed.

## Supplementary figures

**Figure S1. Safety-by-Design weight sensitivity.** **A**, Distribution of Spearman correlations between each alternative formulation ranking and the original ranking across 3,876 component-weight combinations; the dashed line indicates the median correlation. **B**, Median and full range of ranks for the 12 highest-ranked original formulations. Point color denotes the frequency of rank 1, and point size denotes the frequency of top-10 membership.

**Figure S2. Balanced inflammatory-risk distributions and ranking robustness.** **A**, Violin plots show the full balanced record-level distributions of delivery, adaptive activation, innate activation, cytokine burden, off-target activation and Safety-by-Design score in the Low-, Intermediate- and High-risk strata; embedded box plots show the median and interquartile range. **B**, Mean Safety-by-Design scores of the ten leading formulations across risk strata; LNP_0462 is highlighted in black.

**Figure S3. BPCI null comparison and donor-level uncertainty.** **A**, Observed BPCI values compared with the central 95% interval and mean of the exact shuffled-label null at each post-vaccination day. **B**, Median and 95% percentile interval from 2,000 paired-donor bootstrap replicates; point size denotes donor count and open red points show the original unpaired estimator.

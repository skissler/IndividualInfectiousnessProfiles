


  # What the OutbreakTrees database offers

  # The https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3001685 (Taube et al.,
  # PLOS Biology 2022) has 382 transmission trees across 16 diseases. Critically, 137 trees have
  # symptom onset dates — which is exactly what you'd filter for. The data is stored as igraph objects
  # in an https://github.com/DrakeLab/taube-transmission-trees, so it's directly loadable in R.

  # Each tree has nodes (individuals) with attributes including symptom onset dates, and edges
  # representing who-infected-whom. From this you can extract exactly what your likelihood needs:
  # clusters of serial intervals grouped by index case.

  # What your estimation needs vs. what the data provides

  # Your psi_identifiability.R likelihood operates on a list of clusters, where each cluster is a
  # vector of serial intervals (symptom onset of secondary case minus symptom onset of index case) from
  #  the same infector. The data you need per cluster:

  # 1. Who infected whom — the tree edges give this directly
  # 2. Symptom onset dates — to compute serial intervals as onset_secondary - onset_index
  # 3. Grouping by infector — all secondaries of the same index case form one cluster

  # OutbreakTrees gives you all three. For each tree, you'd extract every node that has ≥2 offspring
  # with known onset dates (since your likelihood requires cluster size ≥2 to be informative about
  # within-cluster correlation), compute the serial intervals, and group them by infector.

  # The key challenge: you also need to assume/estimate the incubation period

  # Your model decomposes the serial interval as S = tau + d_j - d_0, where d_0 and d_j are incubation
  # periods (detection delays). The likelihood in psi_identifiability.R requires known a_obs, b_obs —
  # the Gamma parameters of the incubation period distribution. For real data, you'd need to:

  # - Fix incubation period parameters from the literature for each pathogen (these are generally
  # well-estimated — e.g., COVID-19 incubation ~Gamma with mean ~5d, measles ~Gamma with mean ~12d)
  # - Or jointly estimate them, which would be harder and probably unnecessary given how
  # well-characterized incubation periods are

  # This is quite doable. The incubation period literature is mature for most pathogens in the
  # database.

  # What the analysis pipeline would look like

  # 1. Data extraction: Download data_tibble_trees.RDS from GitHub, filter to trees with symptom onset
  # data, extract infector-infectee pairs, compute serial intervals, group by infector into clusters
  # 2. Parameter setup: For each pathogen, pull incubation period Gamma parameters from published
  # estimates
  # 3. Psi estimation: For each pathogen, run compute_psi_posterior() on the real clusters — this gives
  #  you the full posterior over psi. The machinery already exists in psi_identifiability.R
  # 4. Sensitivity analysis: Vary the assumed incubation period parameters to check robustness

  # Practical concerns

  # - Sample size: Your power analysis suggests ~50 index cases suffices for measles and SARS-CoV-2.
  # The question is whether OutbreakTrees has enough clusters per pathogen with onset dates. With 137
  # trees having onset data, and some diseases having many trees (COVID-19 dominates), you may have
  # enough for a few pathogens but not all 16.
  # - Cluster size distribution: Your likelihood gains power from clusters with ≥2 secondaries (sibling
  #  serial intervals reveal within-cluster correlation). Trees dominated by singletons won't help
  # much. The database reportedly has considerable superspreading heterogeneity, which is actually good
  #  for you — superspreaders produce large clusters.
  # - Heterogeneous R0: Your simulation assumes Poisson offspring with fixed R0. Real data will have
  # overdispersed offspring counts. This doesn't break the likelihood (which conditions on observed
  # cluster sizes), but it's worth noting.
  # - Pooling across outbreaks: Different trees for the same pathogen may come from different
  # settings/time periods. You'd need to decide whether to pool them (more data, assumes psi is
  # pathogen-intrinsic) or analyze separately (less power per tree).

  # Other data sources worth considering

  # - Singapore COVID-19 clusters (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/) —
  # well-documented infector-infectee pairs with onset dates, used extensively for generation interval
  # estimation
  # - COVID-19 cluster data from Tianjin, China — same study, second dataset
  # - The eLife COVID-19 infectiousness profile dataset
  # (https://elifesciences.org/articles/79134/figures) — has transmission pair data with onset dates
  # - MERS-CoV hospital clusters — several well-documented nosocomial outbreaks with detailed onset
  # timelines
  # - Mpox 2022 transmission pairs — https://wwwnc.cdc.gov/eid/article/30/6/23-1095_article with serial
  #  interval estimates

  # The OutbreakTrees database is probably the best single starting point because it's already
  # structured as transmission trees across multiple pathogens, and it's directly downloadable as R
  # objects.

  # Bottom line

  # This is quite feasible. The main work is data wrangling (extracting clusters from igraph tree
  # objects), not new methodology — your existing compute_psi_posterior() function is the estimation
  # engine. The biggest uncertainty is whether there are enough large clusters per pathogen in the
  # onset-dated subset of OutbreakTrees.






  # This is very informative. Here's what we're looking at:

  # 11 diseases have ≥5 multi-offspring clusters — enough to potentially estimate psi. The strongest
  # candidates:
  # Disease: MERS
  # Multi-clusters: 23
  # Total SI: 265
  # Notes: Best dataset. One huge cluster (84 offspring from Korean superspreader)
  # ────────────────────────────────────────
  # Disease: COVID-19
  # Multi-clusters: 31
  # Total SI: 191
  # Notes: Good, but pooled across 30 trees from different settings
  # ────────────────────────────────────────
  # Disease: Pneumonic plague
  # Multi-clusters: 20
  # Total SI: 77
  # Notes: Good cluster size distribution
  # ────────────────────────────────────────
  # Disease: Hepatitis A
  # Multi-clusters: 18
  # Total SI: 72
  # Notes: Long serial intervals (~25d mean), many trees
  # ────────────────────────────────────────
  # Disease: Measles
  # Multi-clusters: 8
  # Total SI: 72
  # Notes: Dominated by one large cluster (35 offspring)
  # ────────────────────────────────────────
  # Disease: Norovirus
  # Multi-clusters: 14
  # Total SI: 59
  # Notes: Short serial intervals (~3d)
  # ────────────────────────────────────────
  # Disease: Nipah
  # Multi-clusters: 5
  # Total SI: 55
  # Notes: Two trees, dominated by two big clusters
  # ────────────────────────────────────────
  # Disease: Ebola
  # Multi-clusters: 11
  # Total SI: 43
  # Notes: Pooled from 9 trees
  # ────────────────────────────────────────
  # Disease: Smallpox
  # Multi-clusters: 5
  # Total SI: 46
  # Notes: Single tree, reasonable sizes
  # ────────────────────────────────────────
  # Disease: H1N1
  # Multi-clusters: 9
  # Total SI: 32
  # Notes: Single outbreak, very short SI (~2d)
  # ────────────────────────────────────────
  # Disease: H1N1/H3N2
  # Multi-clusters: 8
  # Total SI: 39
  # Notes: Single outbreak, short SI (~2d)
  # A few things to note:

  # - MERS is the standout — 265 serial intervals with 23 multi-offspring clusters. The Korean 2015
  # outbreak alone has a superspreader with 84 secondary cases, which is a huge cluster for detecting
  # within-cluster correlation.
  # - The data structure maps directly to your likelihood. Each infector with ≥2 offspring gives a
  # cluster of sibling serial intervals, exactly what compute_psi_posterior() needs.
  # - Pooling across trees is a judgment call. For COVID-19, you're pooling across outbreaks in
  # different countries/settings. This is defensible if you argue psi is intrinsic to the pathogen
  # biology (determined by the generation interval distribution shape), but worth flagging.
  # - Incubation period parameters are well-characterized in the literature for most of these
  # (COVID-19, MERS, measles, Ebola, smallpox). That's needed as input (a_obs, b_obs).
  # - The measles data is sparse (only 4 trees, 8 multi-clusters) but you already have measles in your
  # simulation. It would be interesting to see if the real-data estimate is consistent with theory.

  # The analysis runs end-to-end and produces both figures. Here's a summary of the results:

  # Results: Estimated psi across 11 diseases
  # ┌──────────────────┬──────────┬──────────────┬────────────────┬──────────────────┐
  # │     Disease      │ psi mode │    95% CI    │ Multi-clusters │ Serial intervals │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Measles          │ 0.00     │ [0.00, 0.14] │ 8              │ 72               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Pneumonic plague │ 0.00     │ [0.00, 0.66] │ 20             │ 77               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ H1N1             │ 0.15     │ [0.01, 0.82] │ 9              │ 32               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Hepatitis A      │ 0.21     │ [0.06, 0.59] │ 18             │ 72               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Nipah virus      │ 0.23     │ [0.01, 0.97] │ 5              │ 55               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ COVID-19         │ 0.26     │ [0.01, 0.95] │ 31             │ 191              │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ H1N1/H3N2        │ 0.28     │ [0.05, 0.90] │ 8              │ 39               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ MERS             │ 0.43     │ [0.35, 0.53] │ 23             │ 265              │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Smallpox         │ 0.43     │ [0.25, 0.78] │ 5              │ 46               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Norovirus        │ 0.81     │ [0.43, 0.98] │ 14             │ 59               │
  # ├──────────────────┼──────────┼──────────────┼────────────────┼──────────────────┤
  # │ Ebola            │ 1.00     │ [0.03, 0.98] │ 11             │ 43               │
  # └──────────────────┴──────────┴──────────────┴────────────────┴──────────────────┘
  # Key observations:

  # - MERS is the tightest estimate (CI width 0.18) — consistent with having the most data (265 SI, 23
  # multi-clusters). The mode at 0.43 suggests moderately punctuated infectiousness.
  # - Measles strongly favors low psi (mode 0, upper CI 0.14) — punctuated infectiousness, consistent
  # with what you'd expect biologically (narrow generation interval relative to incubation period).
  # - Norovirus favors high psi (mode 0.81) — more diffuse infectiousness, which makes sense for a
  # pathogen with prolonged viral shedding relative to its very short incubation/generation interval.
  # - Wide credible intervals for diseases with fewer multi-offspring clusters (Nipah, COVID-19, Ebola,
  #  H1N1) — these are underpowered, consistent with your power analysis showing ~50+ clusters needed.

  # Important caveats to flag:

  # 1. The GI parameters are moment-matched from observed SI, which is circular in a sense — ideally
  # you'd jointly estimate (alpha, beta, psi) together. For COVID-19, Nipah, and Ebola, the GI variance
  #  had to be floored because the incubation variance explained nearly all SI variance, making the GI
  # estimate fragile.
  # 2. Pooling across trees from different outbreaks/settings assumes psi is intrinsic to the pathogen.
  # 3. The analysis currently treats GI parameters as known rather than jointly estimating them with
  # psi.

  # The script is at code/psi_empirical.R, figures at figures/psi_empirical_posteriors.{pdf,png} and
  # figures/psi_empirical_summary.{pdf,png}, and tabular results at output/psi_empirical_results.csv.
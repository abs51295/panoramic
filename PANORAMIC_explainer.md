# PANORAMIC: What It Does, Why, and the Problem It Solves

## The Problem

Modern spatial omics platforms — Xenium, Visium, MERFISH, CosMx — produce maps of where cells physically sit within tissue. A natural question is: **do two cell types tend to cluster near each other (colocalize), or do they repel?** Such spatial relationships can reveal functional interactions, immune niches, or tumor microenvironments.

The challenge arises the moment you move from a single tissue section to a real clinical study with many patients:

1. **Within-sample uncertainty.** Each tissue section has a finite number of cells. Estimating colocalization from 50 cells per type is inherently noisy.
2. **Between-sample heterogeneity.** Two patients diagnosed with the same disease may have biologically different tumors, so the true colocalization effect can vary across samples — this is expected, not a nuisance.
3. **Group comparison.** The scientific question is usually *differential*: is T-cell / macrophage colocalization stronger in responders versus non-responders? Answering this correctly requires pooling evidence across many samples while accounting for both sources of variability above.

Naive approaches fail here: computing colocalization in each sample separately and then averaging ignores within-sample uncertainty; treating all cells as one big pool ignores sample identity and inflates false positives.

---

## What PANORAMIC Does

**PANORAMIC** (**P**ooled **AN**alysis **O**f Va**R**iance-**A**ware **M**odeling and **I**nference of **C**olocalization) is an R/Bioconductor package implementing a principled three-stage pipeline:

```
Per-sample spatial statistics  →  Random-effects meta-analysis  →  Group comparison
       (Stage 1)                          (Stage 2)                     (Stage 3)
```

It takes a set of `SpatialExperiment` objects (one per patient/sample), each annotated with cell-type labels and XY coordinates, and outputs a ranked list of cell-type pairs that are **differentially colocalized** between two groups, along with effect sizes, standard errors, heterogeneity estimates, and FDR-corrected p-values.

---

## Statistical Framework

### Stage 1: Spatial Statistics with Bootstrap Variance (per sample)

For each pair of cell types (ct1, ct2) and each sample, PANORAMIC estimates the **cross-type L-function** (Ripley's Lcross):

$$L_{12}(r) = \sqrt{\frac{K_{12}(r)}{\pi}}$$

where $K_{12}(r)$ counts the expected number of type-2 cells within distance $r$ of a type-1 cell, normalized by the overall type-2 density. Under **complete spatial randomness (CSR)** — i.e., no association — $L_{12}(r) = r$. The centered version $\hat{L}_{12}(r) - r$ is positive when the two types cluster together and negative when they repel.

**Why the L-function?** It is variance-stabilized relative to the raw K-function, making pooling across scales more straightforward. For self-pairs (ct1 = ct2), the univariate L-function measures within-type clustering.

**Bootstrap variance:** PANORAMIC uses the **Loh bootstrap** (`spatstat.explore::lohboot()`) — a spatial resampling method that respects the point-process structure — to obtain a per-radius, per-sample **variance estimate** $v_i(r)$ alongside the point estimate $y_i(r)$. This variance captures within-sample uncertainty from finite cell counts and spatial configuration.

Values are interpolated onto a common radius grid $r_1, \ldots, r_R$ so that all samples share the same evaluation points.

### Stage 2: Random-Effects Meta-Analysis (pooling across samples)

At each (ct1, ct2, radius) feature, the $n$ sample-level estimates $y_i$ with variances $v_i$ are combined using a **univariate random-effects model** (via `metafor::rma.uni()`):

$$y_i = \mu + u_i + \varepsilon_i, \quad u_i \sim \mathcal{N}(0, \tau^2), \quad \varepsilon_i \sim \mathcal{N}(0, v_i)$$

- $\mu$: the true pooled colocalization at this radius
- $\tau^2$: **between-sample heterogeneity** — how much the true effect varies across patients
- $v_i$: within-sample variance from Stage 1

This model explicitly separates the two sources of variability. The pooled estimate $\hat{\mu}$ is a precision-weighted average that down-weights noisy samples; $\tau^2$ tells you whether the effect is consistent across patients.

**When `group_col` is specified**, separate random-effects models are fitted within each group (e.g., control and case), producing group-specific pooled estimates $\hat{\mu}_{\text{ctrl}}$ and $\hat{\mu}_{\text{case}}$ with their own standard errors and heterogeneity.

**Tau² estimators available:** `"SJ"` (Sidik-Jonkman, robust with small k), `"REML"` (restricted maximum likelihood), `"DL"` (DerSimonian-Laird, classical).

**Output per feature:** $\hat{\mu}$, $\widehat{SE}(\hat{\mu})$, $\hat{\tau}^2$, Cochran's Q, $I^2$, $k$ (number of samples contributing), $p_\mu$.

### Stage 3: Differential Testing (group comparison)

Given group-specific pooled estimates from Stage 2, PANORAMIC computes the **differential colocalization effect**:

$$\hat{\beta}_{\text{diff}} = \hat{\mu}_{\text{case}} - \hat{\mu}_{\text{ctrl}}$$

Under the assumption that group-specific estimates are independent (valid when groups contain different patients):

$$\widehat{SE}(\hat{\beta}_{\text{diff}}) = \sqrt{\widehat{SE}(\hat{\mu}_{\text{ctrl}})^2 + \widehat{SE}(\hat{\mu}_{\text{case}})^2}$$

A **two-sided z-test** yields a p-value, and Benjamini–Hochberg correction controls the false discovery rate (FDR) across all cell-type pairs and radii tested.

**Positive $\hat{\beta}_{\text{diff}}$** means the pair is more colocalized in the case group; **negative** means less colocalized (or more repelled) in cases.

---

## Workflow Step-by-Step

```r
library(panoramic)
library(BiocParallel)

# 1. Prepare: harmonize cell types, build spatial windows, cache spatstat objects
prep_list <- panoramic_prepare(
  spe_list   = my_spe_list,       # Named list of SpatialExperiment objects
  design     = my_design_table,   # data.frame: columns "sample" and "group"
  cell_type  = "cell_type",       # colData column with cell-type labels
  min_cells  = 10,                # Drop rare types (< 10 cells per sample)
  window     = "concave",         # Concave hull for irregular tissue shapes
  concavity  = 50,                # Hull detail (lower = more detail)
  BPPARAM    = MulticoreParam(4)
)

# 2. Compute spatial statistics (Lcross with Loh bootstrap)
se_stats <- panoramic_spatialstats(
  prep_list = prep_list,
  r         = seq(10, 200, by = 10),   # Radius grid in micrometers
  stat      = "Lcross",
  nboot     = 200,
  BPPARAM   = MulticoreParam(4)
)
# Returns SummarizedExperiment:
#   assay "yi": [n_features × n_samples] centered L estimates
#   assay "vi": [n_features × n_samples] bootstrap variances
#   rowData: ct1, ct2, radius_um, stat

# 3. Random-effects meta-analysis per group
se_meta <- panoramic_meta(
  se        = se_stats,
  tau2      = "SJ",          # Sidik-Jonkman tau² estimator
  group_col = "group",       # Fit separate models per group
  BPPARAM   = MulticoreParam(4)
)
# Adds to rowData: {group}_mu_hat, {group}_se_mu, {group}_tau2, {group}_I2, etc.

# 4. Differential comparison between groups
se_diff <- panoramic_compare_groups(
  se     = se_meta,
  group1 = "control",   # Reference group (make.names() applied internally)
  group2 = "case"
)
# Adds to rowData: beta_diff, se_diff, z_diff, p_diff, fdr_diff

# 5. Visualize
plot_volcano(se_diff, fdr_threshold = 0.05, effect_threshold = 1)
plot_forest(se_diff, ct1 = "T_cell", ct2 = "Macrophage", radius = 50)
net <- create_spatial_network(se_diff, fdr_threshold = 0.05)
plot_spatial_network(net)
```

### Convenience wrapper

```r
# Steps 1 + 2 combined:
se_stats <- panoramic(
  spe_list  = my_spe_list,
  design    = my_design_table,
  r         = seq(10, 200, by = 10),
  ...
)
```

---

## Key Data Structures

| Object | Class | Key contents |
|--------|-------|-------------|
| `spe_list[[i]]` | `SpatialExperiment` | Gene expression, cell-type labels in `colData`, XY in `spatialCoords`; spatstat cache in `metadata(spe)$panoramic` |
| `se_stats` | `SummarizedExperiment` | assays `yi` + `vi`; `rowData` has `ct1`, `ct2`, `radius_um`, `stat`; `colData` has `sample`, `group` |
| `se_meta` | `SummarizedExperiment` | Same as above + group-wise meta-analysis columns in `rowData` |
| `se_diff` | `SummarizedExperiment` | Same as above + `beta_diff`, `fdr_diff`, etc. in `rowData` |

---

## Interpreting Results

### Volcano plot
- **X-axis:** Effect size $\hat{\beta}_{\text{diff}}$ (case − control) in centered-L units
- **Y-axis:** $-\log_{10}(p)$
- Points in the upper-right quadrant: pairs **more colocalized** in cases
- Points in the upper-left quadrant: pairs **less colocalized** (or more repelled) in cases
- Color encodes combined FDR + effect-size significance tiers

### Forest plot
Shows individual sample estimates (dots with 95% CI) and the pooled random-effects estimate (diamond) for a single pair at a chosen radius. The width of the diamond reflects precision of the pooled estimate. $\tau^2$ and $I^2$ quantify between-sample heterogeneity.

### Network plot
Cell types are nodes; significant differentially colocalized pairs are edges. Edge width = |z-score|; edge opacity = 1 − FDR. Leiden clustering groups cell types that share similar colocalization patterns. Node size can encode degree, betweenness centrality, or strength.

### Key output columns in `rowData(se_diff)`

| Column | Meaning |
|--------|---------|
| `beta_diff` | Effect size: pooled Lcross(case) − pooled Lcross(control) at this radius |
| `se_diff` | Standard error of `beta_diff` |
| `z_diff` | Z-statistic |
| `p_diff` | Two-sided p-value |
| `fdr_diff` | BH-adjusted p-value across all features |
| `mu_control` / `mu_case` | Group-specific pooled L estimates |
| `tau2_control` / `tau2_case` | Between-patient heterogeneity per group |
| `radius_um` | Radius at which the statistic was evaluated (µm) |

### Heterogeneity interpretation ($\tau^2$ and $I^2$)
- $\tau^2 \approx 0$: the true colocalization is consistent across patients → high confidence in pooled estimate
- $I^2 > 50\%$: substantial heterogeneity → individual patient context matters; the pooled estimate is an average of genuinely variable effects

---

## Why Not a Simpler Approach?

| Approach | Problem |
|----------|---------|
| Pool all cells together | Pseudo-replication: treats one patient with 10,000 cells as 100× more informative than one with 100 cells; inflates N |
| Average per-sample statistics | Ignores within-sample uncertainty; noisy samples contribute equally to precise ones |
| Fixed-effects meta-analysis | Assumes all patients share *exactly* the same true effect; biologically unrealistic |
| Permutation test per sample | Doesn't provide a framework for combining evidence or comparing groups across samples |

PANORAMIC's random-effects meta-analysis is the gold standard from clinical trial literature adapted to spatial biology.

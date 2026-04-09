#!/usr/bin/env Rscript
# =============================================================================
# PANORAMIC differential colocalization analysis for Xenium spatial data
#
# INPUT:  An AnnData (.h5ad) file containing Xenium cells with:
#           obs columns : cell_type, sample_id, group
#           obsm key    : "spatial" (N x 2 matrix of XY coordinates in µm)
#                         Accepted alternatives: "X_spatial", "spatial_xenium"
#                         or obs columns: "x_centroid", "y_centroid"
#
# OUTPUT: results/differential_results.csv
#         results/heatmap_effect_size.pdf  (+ .png)
#         results/radial_profiles_top_pairs.pdf  (+ .png)
#
# USAGE:
#   Rscript xenium_differential_analysis.R \
#     --input      /path/to/data.h5ad \
#     --control    control \
#     --case       case \
#     --outdir     results \
#     --cores      4
#
# DEPENDENCIES (install if missing):
#   BiocManager::install(c("zellkonverter","SpatialExperiment","panoramic",
#                          "BiocParallel","S4Vectors","SummarizedExperiment"))
#   install.packages(c("ggplot2","dplyr","scales","optparse","ggrepel"))
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ---------------------------------------------------------------------------
# 1. Command-line arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--input",       type = "character", default = NULL,
              help = "Path to AnnData .h5ad file [required]"),
  make_option("--control",     type = "character", default = "control",
              help = "Group label for the control/reference group [default: 'control']"),
  make_option("--case",        type = "character", default = "case",
              help = "Group label for the case/comparison group [default: 'case']"),
  make_option("--cell_type_col", type = "character", default = "cell_type",
              help = "obs column containing cell type labels [default: 'cell_type']"),
  make_option("--sample_col",  type = "character", default = "sample_id",
              help = "obs column containing sample IDs [default: 'sample_id']"),
  make_option("--group_col",   type = "character", default = "group",
              help = "obs column containing group labels [default: 'group']"),
  make_option("--spatial_key", type = "character", default = "auto",
              help = "obsm key for spatial coordinates [default: auto-detect]"),
  make_option("--min_cells",   type = "integer",   default = 10L,
              help = "Min cells per type per sample [default: 10]"),
  make_option("--r_min",       type = "double",    default = 10.0,
              help = "Minimum radius in µm [default: 10]"),
  make_option("--r_max",       type = "double",    default = 200.0,
              help = "Maximum radius in µm [default: 200]"),
  make_option("--r_step",      type = "double",    default = 10.0,
              help = "Radius step size in µm [default: 10]"),
  make_option("--nsim",        type = "integer",   default = 199L,
              help = "Bootstrap simulations for variance estimation [default: 199]"),
  make_option("--fdr_threshold", type = "double",  default = 0.05,
              help = "FDR significance threshold for plots [default: 0.05]"),
  make_option("--top_pairs",   type = "integer",   default = 10L,
              help = "Number of top pairs shown in radial profile plot [default: 10]"),
  make_option("--outdir",      type = "character", default = "results",
              help = "Output directory [default: 'results']"),
  make_option("--cores",       type = "integer",   default = 1L,
              help = "Number of CPU cores for parallel computation [default: 1]"),
  make_option("--seed",        type = "integer",   default = 42L,
              help = "Random seed [default: 42]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("--input is required. Provide a path to an .h5ad file.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 2. Load required packages
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(BiocParallel)
  library(panoramic)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

bpparam <- if (opt$cores > 1L) {
  BiocParallel::MulticoreParam(opt$cores, RNGseed = opt$seed)
} else {
  BiocParallel::SerialParam()
}

message("=== PANORAMIC Xenium Differential Colocalization ===")
message("Input : ", opt$input)
message("Groups: ", opt$control, " vs ", opt$case)
message("Radii : ", opt$r_min, " – ", opt$r_max, " µm (step ", opt$r_step, " µm)")
message("Cores : ", opt$cores)

# ---------------------------------------------------------------------------
# 3. Load AnnData and convert to list of SpatialExperiment objects
# ---------------------------------------------------------------------------

load_anndata_as_spe_list <- function(
    h5ad_path,
    cell_type_col = "cell_type",
    sample_col    = "sample_id",
    group_col     = "group",
    spatial_key   = "auto"
) {
  message("Loading AnnData from: ", h5ad_path)
  sce <- zellkonverter::readH5AD(h5ad_path, use_hdf5 = FALSE, verbose = FALSE)

  # ----- Resolve spatial coordinates -----
  find_coords <- function(sce, key) {
    # 1. Try obsm keys stored as reducedDims
    rdn <- SingleCellExperiment::reducedDimNames(sce)

    candidates <- if (key == "auto") {
      c("spatial", "X_spatial", "spatial_xenium", "Spatial")
    } else {
      key
    }

    coords <- NULL
    for (k in candidates) {
      if (k %in% rdn) {
        mat <- SingleCellExperiment::reducedDim(sce, k)
        if (ncol(mat) >= 2) {
          coords <- mat[, 1:2, drop = FALSE]
          colnames(coords) <- c("x", "y")
          message("  Spatial coordinates from reducedDim('", k, "')")
          break
        }
      }
    }

    # 2. Fallback: look for x_centroid / y_centroid in colData (10x Xenium format)
    if (is.null(coords)) {
      cd <- as.data.frame(SummarizedExperiment::colData(sce))
      xy_cols <- intersect(c("x_centroid", "y_centroid"), colnames(cd))
      if (length(xy_cols) == 2) {
        coords <- as.matrix(cd[, xy_cols])
        colnames(coords) <- c("x", "y")
        message("  Spatial coordinates from colData columns: x_centroid, y_centroid")
      }
    }

    if (is.null(coords)) {
      stop(
        "Cannot find spatial coordinates. Checked reducedDims: ",
        paste(rdn, collapse = ", "),
        ". Also checked colData for x_centroid/y_centroid. ",
        "Pass --spatial_key to specify the obsm key explicitly."
      )
    }
    storage.mode(coords) <- "double"
    coords
  }

  coords <- find_coords(sce, spatial_key)

  # ----- Validate required colData columns -----
  cd  <- as.data.frame(SummarizedExperiment::colData(sce))
  for (col in c(cell_type_col, sample_col, group_col)) {
    if (!col %in% colnames(cd)) {
      stop("Required obs column '", col, "' not found in the AnnData. ",
           "Available columns: ", paste(colnames(cd), collapse = ", "))
    }
  }

  cell_types <- as.character(cd[[cell_type_col]])
  sample_ids <- as.character(cd[[sample_col]])
  groups     <- as.character(cd[[group_col]])

  message("  Cells        : ", nrow(cd))
  message("  Cell types   : ", paste(sort(unique(cell_types)), collapse = ", "))
  message("  Samples      : ", paste(sort(unique(sample_ids)), collapse = ", "))
  message("  Groups       : ", paste(sort(unique(groups)), collapse = ", "))

  # ----- Split into per-sample SpatialExperiment objects -----
  sample_levels <- unique(sample_ids)
  spe_list <- lapply(sample_levels, function(sid) {
    idx <- which(sample_ids == sid)

    # Build a minimal counts matrix (1 x n_cells) as a placeholder
    # (PANORAMIC uses cell coordinates + types, not expression)
    cnt <- matrix(1L, nrow = 1L, ncol = length(idx),
                  dimnames = list("placeholder", paste0(sid, "_", seq_along(idx))))

    SpatialExperiment::SpatialExperiment(
      assays        = list(counts = cnt),
      colData       = S4Vectors::DataFrame(
        cell_type = cell_types[idx],
        sample_id = sid,
        group     = groups[idx][1]   # group is constant within a sample
      ),
      spatialCoords = coords[idx, , drop = FALSE]
    )
  })
  names(spe_list) <- sample_levels

  # ----- Build design table -----
  design <- data.frame(
    sample = sample_levels,
    group  = vapply(sample_levels,
                    function(sid) unique(groups[sample_ids == sid])[1],
                    character(1)),
    stringsAsFactors = FALSE
  )

  list(spe_list = spe_list, design = design)
}

data_obj <- load_anndata_as_spe_list(
  h5ad_path     = opt$input,
  cell_type_col = opt$cell_type_col,
  sample_col    = opt$sample_col,
  group_col     = opt$group_col,
  spatial_key   = opt$spatial_key
)
spe_list <- data_obj$spe_list
design   <- data_obj$design

# Verify both groups are present
groups_found <- unique(design$group)
for (g in c(opt$control, opt$case)) {
  if (!g %in% groups_found) {
    stop("Group '", g, "' not found in the data. Found: ",
         paste(groups_found, collapse = ", "))
  }
}

# ---------------------------------------------------------------------------
# 4. PANORAMIC pipeline – Xenium-optimized parameters
# ---------------------------------------------------------------------------

# Radius grid (10–200 µm captures 1–20 cell-diameter neighborhoods)
radii_um <- seq(opt$r_min, opt$r_max, by = opt$r_step)

message("\n--- Stage 1: Prepare spatial windows ---")
prep_list <- panoramic_prepare(
  spe_list  = spe_list,
  design    = design,
  cell_type = "cell_type",
  min_cells = opt$min_cells,    # Drop types with < min_cells per sample
  window    = "concave",        # Concave hull for irregular tissue sections
  concavity = 50,               # Moderate detail; increase for tighter hulls
  BPPARAM   = bpparam
)

message("\n--- Stage 2: Spatial statistics (Lcross + Loh bootstrap) ---")
se_stats <- panoramic_spatialstats(
  prep      = prep_list,
  pairs     = "auto",           # All observed cell-type pairs
  radii_um  = radii_um,
  stat      = "Lcross",         # Cross-type L-function (Ripley's)
  nsim      = opt$nsim,         # Bootstrap samples for variance
  correction = "translate",     # Edge correction (robust for tissue sections)
  seed      = opt$seed,
  BPPARAM   = bpparam,
  verbose   = FALSE
)
message("  Features: ", nrow(se_stats),
        " (", length(unique(rowData(se_stats)$ct1)), " x ",
        length(unique(rowData(se_stats)$ct2)), " pairs x ",
        length(radii_um), " radii)")

message("\n--- Stage 3: Random-effects meta-analysis per group ---")
se_meta <- panoramic_meta(
  se        = se_stats,
  tau2      = "SJ",             # Sidik-Jonkman: robust when k is small
  group_col = "group",
  BPPARAM   = bpparam
)

message("\n--- Stage 4: Differential comparison (", opt$case, " vs ", opt$control, ") ---")
# make.names() is applied internally by panoramic_meta to group labels
ctrl_name <- make.names(opt$control)
case_name <- make.names(opt$case)

se_diff <- panoramic_compare_groups(
  se     = se_meta,
  group1 = ctrl_name,
  group2 = case_name
)

# ---------------------------------------------------------------------------
# 5. Save results table
# ---------------------------------------------------------------------------
rd <- as.data.frame(SummarizedExperiment::rowData(se_diff))
results_path <- file.path(opt$outdir, "differential_results.csv")
write.csv(rd, results_path, row.names = FALSE)
message("\nResults saved to: ", results_path)

n_sig <- sum(!is.na(rd$fdr_diff) & rd$fdr_diff < opt$fdr_threshold, na.rm = TRUE)
message("Significant pairs (FDR < ", opt$fdr_threshold, "): ", n_sig,
        " features out of ", nrow(rd))

# ---------------------------------------------------------------------------
# 6. Heatmap: effect size and significance across cell-type pairs
# ---------------------------------------------------------------------------
# For each (ct1, ct2) pair, summarize across radii by selecting the radius
# with the most significant (smallest p_diff) result. This gives a single
# representative effect size per pair that captures the most informative scale.

make_heatmap <- function(rd, fdr_threshold, title_suffix = "") {

  # Per-pair summary: pick the radius of peak significance
  pair_summary <- rd %>%
    filter(!is.na(beta_diff) & !is.na(fdr_diff)) %>%
    group_by(ct1, ct2) %>%
    slice_min(order_by = fdr_diff, n = 1, with_ties = FALSE) %>%
    ungroup()

  if (nrow(pair_summary) == 0) {
    message("No valid differential results to plot in heatmap.")
    return(NULL)
  }

  # Significance label
  pair_summary <- pair_summary %>%
    mutate(
      sig_label = case_when(
        fdr_diff < 0.001 ~ "***",
        fdr_diff < 0.01  ~ "** ",
        fdr_diff < 0.05  ~ "*  ",
        TRUE             ~ ""
      )
    )

  # Symmetric axis ordering: sort all cell types alphabetically
  all_cts <- sort(union(pair_summary$ct1, pair_summary$ct2))

  pair_summary <- pair_summary %>%
    mutate(
      ct1 = factor(ct1, levels = all_cts),
      ct2 = factor(ct2, levels = rev(all_cts))
    )

  # Color scale limits: symmetric around 0
  max_abs <- max(abs(pair_summary$beta_diff), na.rm = TRUE)
  lim <- ceiling(max_abs * 10) / 10

  p <- ggplot(pair_summary, aes(x = ct2, y = ct1)) +
    geom_tile(aes(fill = beta_diff), color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = sig_label),
      size = 3.5, color = "black", vjust = 0.5, family = "mono"
    ) +
    scale_fill_gradientn(
      colours  = c("#053061", "#2166ac", "#92c5de", "#f7f7f7",
                   "#f4a582", "#d6604d", "#67001f"),
      limits   = c(-lim, lim),
      oob      = scales::squish,
      name     = expression(hat(beta)[diff]),
      guide    = guide_colorbar(barheight = 12, title.position = "top",
                                title.hjust = 0.5)
    ) +
    labs(
      x        = "ct2 (neighbor)",
      y        = "ct1 (reference)",
      title    = paste0("Differential spatial colocalization: ",
                        title_suffix),
      subtitle = paste0(
        "Fill = effect size (case \u2212 control); ",
        "***\u202fFDR<0.001, **\u202fFDR<0.01, *\u202fFDR<0.05\n",
        "Radius shown: peak significance per pair"
      ),
      caption  = paste0(
        "n = ", length(unique(pair_summary$ct1)), " cell types; ",
        "FDR threshold = ", fdr_threshold
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y      = element_text(size = 9),
      axis.title       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 9, color = "gray40"),
      legend.position  = "right",
      panel.grid       = element_blank()
    )

  p
}

title_sfx <- paste0(opt$case, " vs ", opt$control)
p_heatmap <- make_heatmap(rd, fdr_threshold = opt$fdr_threshold,
                          title_suffix = title_sfx)

if (!is.null(p_heatmap)) {
  n_cts <- length(union(unique(rd$ct1), unique(rd$ct2)))
  hw    <- max(5, n_cts * 0.55)   # scale figure with number of cell types

  heatmap_base <- file.path(opt$outdir, "heatmap_effect_size")
  ggsave(paste0(heatmap_base, ".pdf"), p_heatmap,
         width = hw + 2, height = hw, device = cairo_pdf)
  ggsave(paste0(heatmap_base, ".png"), p_heatmap,
         width = hw + 2, height = hw, dpi = 150)
  message("Heatmap saved to: ", heatmap_base, ".pdf / .png")
}

# ---------------------------------------------------------------------------
# 7. Radial profiles: effect size vs radius for top significant pairs
# ---------------------------------------------------------------------------
make_radial_profiles <- function(rd, top_n = 10, fdr_threshold = 0.05,
                                  title_suffix = "") {

  # Rank pairs by minimum FDR across all radii
  pair_ranks <- rd %>%
    filter(!is.na(fdr_diff)) %>%
    group_by(ct1, ct2) %>%
    summarise(min_fdr = min(fdr_diff, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_fdr) %>%
    slice_head(n = top_n)

  if (nrow(pair_ranks) == 0) {
    message("No pairs with valid FDR for radial profile plot.")
    return(NULL)
  }

  # Filter rowData to these pairs
  plot_rd <- rd %>%
    inner_join(pair_ranks %>% select(ct1, ct2), by = c("ct1", "ct2")) %>%
    filter(!is.na(beta_diff)) %>%
    mutate(
      pair_label = paste0(ct1, "\u2192", ct2),
      significant = !is.na(fdr_diff) & fdr_diff < fdr_threshold
    )

  # Order by min FDR
  ordered_labels <- paste0(pair_ranks$ct1, "\u2192", pair_ranks$ct2)
  plot_rd$pair_label <- factor(plot_rd$pair_label, levels = ordered_labels)

  # Color points by significance
  p <- ggplot(plot_rd, aes(x = radius_um, y = beta_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_ribbon(
      aes(ymin = beta_diff - 1.96 * se_diff,
          ymax = beta_diff + 1.96 * se_diff),
      fill = "#4393c3", alpha = 0.15
    ) +
    geom_line(color = "#2166ac", linewidth = 0.8) +
    geom_point(
      aes(color = significant, shape = significant),
      size = 2.2
    ) +
    scale_color_manual(
      values = c("TRUE" = "#d62728", "FALSE" = "gray60"),
      labels = c("TRUE" = paste0("FDR < ", fdr_threshold), "FALSE" = "ns"),
      name   = "Significance"
    ) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = paste0("FDR < ", fdr_threshold), "FALSE" = "ns"),
      name   = "Significance"
    ) +
    facet_wrap(~ pair_label, ncol = 2, scales = "free_y") +
    labs(
      x        = "Radius (µm)",
      y        = expression(hat(beta)[diff] ~ "(case \u2212 control)"),
      title    = paste0("Radial profiles of top ", nrow(pair_ranks),
                        " pairs: ", title_suffix),
      subtitle = paste0("Shaded band: 95% CI; red points: FDR < ",
                        fdr_threshold),
      caption  = "Pairs ranked by minimum FDR across radii"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text    = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "#e8f4f8"),
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "bottom"
    )

  p
}

p_radial <- make_radial_profiles(
  rd,
  top_n         = opt$top_pairs,
  fdr_threshold = opt$fdr_threshold,
  title_suffix  = title_sfx
)

if (!is.null(p_radial)) {
  n_top      <- min(opt$top_pairs, nrow(rd))
  n_cols     <- 2L
  n_rows_fig <- ceiling(n_top / n_cols)
  pw <- 10
  ph <- max(4, n_rows_fig * 3.5)

  radial_base <- file.path(opt$outdir, "radial_profiles_top_pairs")
  ggsave(paste0(radial_base, ".pdf"), p_radial,
         width = pw, height = ph, device = cairo_pdf)
  ggsave(paste0(radial_base, ".png"), p_radial,
         width = pw, height = ph, dpi = 150)
  message("Radial profiles saved to: ", radial_base, ".pdf / .png")
}

# ---------------------------------------------------------------------------
# 8. Standard PANORAMIC volcano plot (bonus)
# ---------------------------------------------------------------------------
if (n_sig > 0) {
  p_vol <- tryCatch(
    panoramic::plot_volcano(
      se_diff,
      fdr_threshold    = opt$fdr_threshold,
      effect_threshold = 1.0,
      label_top        = 15
    ),
    error = function(e) {
      message("Volcano plot skipped: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(p_vol)) {
    volcano_base <- file.path(opt$outdir, "volcano_plot")
    ggsave(paste0(volcano_base, ".pdf"), p_vol,
           width = 9, height = 7, device = cairo_pdf)
    ggsave(paste0(volcano_base, ".png"), p_vol,
           width = 9, height = 7, dpi = 150)
    message("Volcano plot saved to: ", volcano_base, ".pdf / .png")
  }
}

message("\n=== Analysis complete. Outputs in: ", opt$outdir, " ===")

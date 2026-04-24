# ============================================================
# Reproducible sPLS-DA workflow with stable keepX selection
# HC / BCN / BCP comparison
# FULLY QUALIFIED WITH package::function()
# ============================================================

# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

suppressWarnings({
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    library(BiocParallel)
  }
})

# ---- [Global Reproducibility Seed] ----
SEED <- 1
set.seed(SEED)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots_P257_09"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Load Expression Data] ----
data <- utils::read.table(
  "08.12.25_NEW_Adjusted_Normalized_sepsis_P257_09.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- [Fix Gene Names] ----
data$Gene_names <- base::ifelse(
  base::is.na(data$Gene_names) | base::trimws(data$Gene_names) == "",
  data$Protein_IDs,
  data$Gene_names
)

data$Gene_names <- base::make.unique(base::as.character(data$Gene_names))

# Use short safe feature IDs for modeling
data$Safe_Feature_ID <- base::paste0("Feature_", base::seq_len(nrow(data)))
rownames(data) <- data$Safe_Feature_ID

# ---- [Define Groups] ----
group_cols <- list(
  HC  = base::grep("^HC_[1-6]$", colnames(data), value = TRUE),
  BCN = base::grep("^BCN_[1-6]$", colnames(data), value = TRUE),
  BCP = base::grep("^BCP_[1-6]$", colnames(data), value = TRUE)
)

selected_groups <- names(group_cols)

group_sizes <- base::sapply(group_cols, length)
print(group_sizes)

if (any(group_sizes == 0)) {
  stop("One or more groups have zero matched columns. Please check sample names.")
}

# ---- [Create X and Y] ----
X_data <- data[, unlist(group_cols), drop = FALSE]
X <- base::as.data.frame(base::t(X_data))

sample_names <- rownames(X)
Y <- dplyr::case_when(
  base::grepl("^HC", sample_names) ~ "HC",
  base::grepl("^BCN", sample_names) ~ "BCN",
  base::grepl("^BCP", sample_names) ~ "BCP",
  TRUE ~ NA_character_
)

valid_idx <- which(!base::is.na(Y))
if (length(valid_idx) == 0) {
  stop("No valid samples matched the defined groups.")
}

X <- X[valid_idx, , drop = FALSE]
Y <- base::factor(Y[valid_idx], levels = selected_groups)

X <- base::as.data.frame(base::lapply(X, base::as.numeric))
rownames(X) <- sample_names[valid_idx]

# Remove columns with all NA or zero variance
nzv <- base::sapply(X, function(z) {
  z <- z[!base::is.na(z)]
  if (length(z) <= 1) return(FALSE)
  stats::sd(z) > 0
})

cat("Variables removed (all NA / zero variance):", sum(!nzv), "\n")
X <- X[, nzv, drop = FALSE]

if (ncol(X) == 0) {
  stop("No valid variables remain after removing all-NA and zero-variance columns.")
}

cat("Samples:", nrow(X), "\n")
cat("Variables after zero-variance filtering:", ncol(X), "\n")
cat("Class counts:\n")
print(base::table(Y))

# ---- [Color Palette] ----
group_colors <- c(
  "HC"  = "#E41A1C",
  "BCN" = "#377EB8",
  "BCP" = "#4DAF4A"
)

# ---- [Custom Theme] ----
custom_theme_splsda <- function(base_size = 25) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size),
      legend.position = "none"
    )
}

# ---- [Tuning Parameters] ----
optimal_ncomp <- 2
list_keepX <- c(25, 50, 100)

bp_param <- NULL
if ("BiocParallel" %in% loadedNamespaces()) {
  bp_param <- BiocParallel::SerialParam()
}

# ---- [Reproducible Tuning sPLS-DA] ----
cat("Tuning sPLS-DA with fixed seed...\n")
set.seed(SEED)

tune_args <- list(
  X = X,
  Y = Y,
  ncomp = optimal_ncomp,
  validation = "Mfold",
  folds = 5,
  dist = "centroids.dist",
  measure = "BER",
  test.keepX = list_keepX,
  nrepeat = 10,
  progressBar = TRUE
)

tune_formals <- names(base::formals(mixOmics::tune.splsda))
if ("seed" %in% tune_formals) {
  tune_args$seed <- SEED
}

if (!base::is.null(bp_param) && "BPPARAM" %in% tune_formals) {
  tune_args$BPPARAM <- bp_param
}

tune_result <- base::do.call(mixOmics::tune.splsda, tune_args)

optimal_keepX <- tune_result$choice.keepX[seq_len(optimal_ncomp)]
print(optimal_keepX)

# ---- [Save Tuning Summary] ----
tuning_summary <- base::data.frame(
  Component = base::paste0("Comp", base::seq_along(optimal_keepX)),
  keepX = base::as.numeric(optimal_keepX)
)

utils::write.table(
  tuning_summary,
  file = file.path(output_dir, "optimal_keepX_summary_P257_09.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---- [keepX Barplot] ----
keepx_df <- base::data.frame(
  Component = base::factor(
    base::paste0("Comp", base::seq_along(optimal_keepX)),
    levels = base::paste0("Comp", base::seq_along(optimal_keepX))
  ),
  keepX = base::as.numeric(optimal_keepX)
)

gg_keepx <- ggplot2::ggplot(keepx_df, ggplot2::aes(x = Component, y = keepX)) +
  ggplot2::geom_col(fill = "#377EB8", alpha = 0.85, width = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = keepX), vjust = -0.5, size = 5) +
  ggplot2::labs(
    title = "Optimal keepX per Component",
    x = "Component",
    y = "Variables Selected"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA)
  )

print(gg_keepx)

ggplot2::ggsave(
  file.path(output_dir, "keepX_per_component_P257_09.png"),
  plot = gg_keepx,
  dpi = 300,
  width = 6,
  height = 5,
  bg = "white"
)

# ---- [Fit Final Model] ----
set.seed(SEED)
splsda_model <- mixOmics::splsda(
  X = X,
  Y = Y,
  ncomp = optimal_ncomp,
  keepX = optimal_keepX
)

# ---- [NEW: Explicit Cross-validated Performance Assessment] ----
cat("\nRunning repeated cross-validated performance assessment...\n")
set.seed(SEED)

perf_args <- list(
  object = splsda_model,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  dist = "centroids.dist",
  progressBar = TRUE
)

perf_formals <- names(base::formals(mixOmics::perf))
if ("seed" %in% perf_formals) {
  perf_args$seed <- SEED
}

if (!base::is.null(bp_param) && "BPPARAM" %in% perf_formals) {
  perf_args$BPPARAM <- bp_param
}

perf_result <- base::do.call(mixOmics::perf, perf_args)

# ---- [NEW: Extract BER summary] ----
ber_by_comp <- NULL

if (!is.null(perf_result$error.rate$BER)) {
  ber_by_comp <- perf_result$error.rate$BER
} else if (!is.null(perf_result$error.rate$centroids.dist)) {
  ber_by_comp <- perf_result$error.rate$centroids.dist
}

if (is.null(ber_by_comp)) {
  warning("Could not extract BER directly from perf_result$error.rate.")
} else {
  ber_summary <- base::data.frame(
    Component = base::paste0("Comp", seq_along(ber_by_comp)),
    BER = base::as.numeric(ber_by_comp),
    stringsAsFactors = FALSE
  )
  
  utils::write.table(
    x = ber_summary,
    file = file.path(output_dir, "BER_by_component_P257_09.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  print(ber_summary)
}

# ---- [NEW: Confusion Matrix from Final Number of Components] ----
conf_mat <- NULL

if (!is.null(perf_result$confusion)) {
  if (!is.null(perf_result$confusion$centroids.dist)) {
    conf_obj <- perf_result$confusion$centroids.dist
    
    if (length(dim(conf_obj)) == 3) {
      conf_mat <- conf_obj[, , optimal_ncomp, drop = FALSE]
      conf_mat <- conf_mat[, , 1]
    } else {
      conf_mat <- conf_obj
    }
  } else {
    conf_mat <- perf_result$confusion
  }
}

if (!is.null(conf_mat)) {
  conf_mat_df <- base::as.data.frame.matrix(conf_mat)
  
  utils::write.table(
    x = cbind(True_Class = rownames(conf_mat_df), conf_mat_df),
    file = file.path(output_dir, "confusion_matrix_P257_09.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  print(conf_mat)
} else {
  warning("Confusion matrix could not be extracted automatically from perf_result.")
}

# ---- [NEW: Optional BER barplot] ----
if (exists("ber_summary")) {
  gg_ber <- ggplot2::ggplot(ber_summary, ggplot2::aes(x = Component, y = BER)) +
    ggplot2::geom_col(fill = "#D95F02", alpha = 0.85, width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = round(BER, 3)), vjust = -0.5, size = 5) +
    ggplot2::labs(
      title = "Cross-validated BER per Component",
      x = "Component",
      y = "Balanced Error Rate"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
  
  print(gg_ber)
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "BER_per_component_P257_09.png"),
    plot = gg_ber,
    dpi = 300,
    width = 6,
    height = 5,
    bg = "white"
  )
}

# ---- [Explained Variance] ----
expl_var <- base::apply(splsda_model$variates$X^2, 2, sum, na.rm = TRUE) /
  sum(splsda_model$variates$X^2, na.rm = TRUE)

# ---- [Robust Ellipse Function] ----
desired_ellipse_level <- 0.95

compute_ellipse <- function(mean_vec, cov_mat, level = 0.95, npoints = 200) {
  cov_mat <- base::as.matrix(cov_mat)
  if (any(!base::is.finite(cov_mat))) return(NULL)
  
  eps <- 1e-8
  cov_mat <- cov_mat + base::diag(eps, nrow(cov_mat))
  
  eig <- base::eigen(cov_mat)
  eig$values[eig$values < eps] <- eps
  
  angles <- base::seq(0, 2 * base::pi, length.out = npoints)
  circle <- base::rbind(base::cos(angles), base::sin(angles))
  radius <- base::sqrt(stats::qchisq(level, df = 2))
  
  transform_mat <- eig$vectors %*% base::diag(base::sqrt(eig$values))
  ellipse <- base::t(radius * transform_mat %*% circle)
  ellipse <- base::sweep(ellipse, 2, mean_vec, "+")
  
  df <- base::as.data.frame(ellipse)
  colnames(df) <- c("comp1", "comp2")
  df
}

# ---- [Plot All Component Pairs] ----
for (i in 1:(optimal_ncomp - 1)) {
  for (j in (i + 1):optimal_ncomp) {
    
    plot_data <- base::data.frame(
      comp1 = splsda_model$variates$X[, i],
      comp2 = splsda_model$variates$X[, j],
      Group = Y,
      Sample = rownames(X),
      stringsAsFactors = FALSE
    )
    
    group_centroids <- plot_data %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        comp1 = mean(comp1, na.rm = TRUE),
        comp2 = mean(comp2, na.rm = TRUE),
        count = dplyr::n(),
        .groups = "drop"
      )
    
    ellipse_list <- base::lapply(base::split(plot_data, plot_data$Group), function(df_group) {
      if (nrow(df_group) < 3) return(NULL)
      
      group_data <- df_group[, c("comp1", "comp2"), drop = FALSE]
      cov_mat <- stats::cov(group_data)
      
      ell <- base::tryCatch(
        compute_ellipse(
          mean_vec = base::colMeans(group_data),
          cov_mat = cov_mat,
          level = desired_ellipse_level
        ),
        error = function(e) NULL
      )
      
      if (base::is.null(ell)) return(NULL)
      ell$Group <- unique(df_group$Group)
      ell
    })
    
    ellipse_data <- dplyr::bind_rows(ellipse_list)
    
    x_lab <- base::paste0("Component ", i, " (", round(expl_var[i] * 100, 1), "%)")
    y_lab <- base::paste0("Component ", j, " (", round(expl_var[j] * 100, 1), "%)")
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = comp1, y = comp2, color = Group, fill = Group)
    ) +
      ggplot2::geom_point(size = 4, alpha = 0.9) +
      {
        if (nrow(ellipse_data) > 0) {
          ggplot2::geom_polygon(
            data = ellipse_data,
            ggplot2::aes(group = Group),
            alpha = 0.2,
            color = NA
          )
        }
      } +
      {
        if (nrow(ellipse_data) > 0) {
          ggplot2::geom_path(
            data = ellipse_data,
            ggplot2::aes(group = Group),
            linewidth = 1
          )
        }
      } +
      ggrepel::geom_text_repel(
        data = group_centroids,
        ggplot2::aes(label = base::paste0(Group, "\n(n=", count, ")")),
        color = "black",
        size = 7,
        fontface = "bold",
        max.overlaps = 100,
        box.padding = 0.6,
        point.padding = 0.6,
        seed = SEED
      ) +
      ggplot2::scale_color_manual(values = group_colors) +
      ggplot2::scale_fill_manual(values = group_colors) +
      ggplot2::labs(
        title = base::paste0("sPLS-DA: P257_09 (Comp ", i, " vs ", j, ")"),
        x = x_lab,
        y = y_lab
      ) +
      custom_theme_splsda()
    
    print(p)
    
    ggplot2::ggsave(
      filename = file.path(output_dir, base::paste0("sPLS-DA_P257_09_Comp", i, "_vs_Comp", j, ".png")),
      plot = p,
      dpi = 300,
      width = 10,
      height = 10,
      bg = "white"
    )
  }
}

# ---- [Session/Version Log for Reproducibility] ----
base::sink(file.path(output_dir, "sessionInfo_P257_09.txt"))
cat("Seed used:", SEED, "\n\n")
print(utils::sessionInfo())
base::sink()

# ---- [Summary Report] ----
cat("\nsPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("Seed Used:", SEED, "\n")
cat("Number of Components Used:", optimal_ncomp, "\n")
cat("keepX per Component:", base::paste(optimal_keepX, collapse = ", "), "\n")
cat("Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")

if (exists("ber_summary")) {
  cat("Cross-validated BER by component:\n")
  print(ber_summary)
}

n_to_report <- min(2, length(expl_var))
cat(
  "Explained Variance (First ", n_to_report, " components): ",
  base::paste0(round(expl_var[1:n_to_report] * 100, 1), collapse = "%, "),
  "%\n",
  sep = ""
)

cat("Output folder:", normalizePath(output_dir), "\n")
cat("\nAll sPLS-DA 2D plots and validation summaries saved successfully.\n")


################################################################################
# STEP 2
# Export top 20 contributing proteins per sPLS-DA component to TSV
################################################################################

library(readr)

top_n <- 20
n_components_to_report <- min(2, ncol(splsda_model$loadings$X))
output_dir_top <- "Top_Contributing_Proteins_P257_09"
dir.create(output_dir_top, showWarnings = FALSE, recursive = TRUE)

if (!exists("data")) {
  stop("Object `data` not found. Please run STEP 1 first.")
}
if (!exists("splsda_model")) {
  stop("Object `splsda_model` not found. Please run STEP 1 first.")
}

data_export <- data

data_export$Gene_names <- base::ifelse(
  base::is.na(data_export$Gene_names) | base::trimws(data_export$Gene_names) == "",
  data_export$Protein_IDs,
  data_export$Gene_names
)

data_export$Gene_names <- base::make.unique(base::as.character(data_export$Gene_names))
data_export$Feature_ID_for_join <- data_export$Safe_Feature_ID

all_data_columns <- colnames(data_export)

top_tables <- base::lapply(base::seq_len(n_components_to_report), function(comp) {
  
  loadings_vec <- splsda_model$loadings$X[, comp]
  loadings_vec <- loadings_vec[!base::is.na(loadings_vec)]
  
  ord <- base::order(base::abs(loadings_vec), decreasing = TRUE)
  top_features <- loadings_vec[ord][base::seq_len(min(top_n, length(loadings_vec)))]
  
  feature_df <- base::data.frame(
    Feature_ID_for_join = names(top_features),
    Loading = base::as.numeric(top_features),
    Abs_Loading = base::abs(base::as.numeric(top_features)),
    Rank = base::seq_along(top_features),
    Component = base::paste0("Comp", comp),
    stringsAsFactors = FALSE
  )
  
  merged_df <- feature_df %>%
    dplyr::left_join(data_export, by = "Feature_ID_for_join") %>%
    dplyr::select(
      Component,
      Rank,
      Gene_names,
      Feature_ID_for_join,
      Loading,
      Abs_Loading,
      dplyr::all_of(all_data_columns)
    )
  
  merged_df
})

top_contributing_all <- dplyr::bind_rows(top_tables)

unmatched_df <- top_contributing_all %>%
  dplyr::filter(base::is.na(Protein_IDs))

if (nrow(unmatched_df) > 0) {
  utils::write.table(
    unmatched_df,
    file = file.path(output_dir_top, "Unmatched_Top_Features.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  cat("Warning:", nrow(unmatched_df), "top features could not be matched back to `data`.\n")
}

utils::write.table(
  top_contributing_all,
  file = file.path(output_dir_top, "Top20_Contributing_Proteins_All_Components.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

for (comp in unique(top_contributing_all$Component)) {
  comp_df <- top_contributing_all %>%
    dplyr::filter(Component == comp)
  
  utils::write.table(
    comp_df,
    file = file.path(output_dir_top, base::paste0("Top20_Contributing_", comp, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

summary_cols <- c(
  "Component",
  "Rank",
  "Gene_names",
  "Feature_ID_for_join",
  "Protein_IDs",
  "Majority_protein_IDs",
  "Protein_names",
  "Loading",
  "Abs_Loading"
)

summary_cols <- summary_cols[summary_cols %in% colnames(top_contributing_all)]

summary_df <- top_contributing_all %>%
  dplyr::select(dplyr::all_of(summary_cols))

utils::write.table(
  summary_df,
  file = file.path(output_dir_top, "Top20_Contributing_Proteins_Summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

component_colors <- c(
  "Comp1" = "#F8766D",
  "Comp2" = "#7CAE00",
  "Comp3" = "#00BFC4",
  "Comp4" = "#C77CFF",
  "Comp5" = "#FF61C3",
  "Comp6" = "#00BA38"
)

plot_df <- top_contributing_all %>%
  dplyr::mutate(
    Label = dplyr::case_when(
      !base::is.na(Gene_names) & base::trimws(Gene_names) != "" ~ Gene_names,
      !base::is.na(Protein_IDs) & base::trimws(Protein_IDs) != "" ~ Protein_IDs,
      TRUE ~ Feature_ID_for_join
    )
  )

for (comp in unique(plot_df$Component)) {
  comp_df <- plot_df %>%
    dplyr::filter(Component == comp) %>%
    dplyr::arrange(Abs_Loading) %>%
    dplyr::mutate(Label = base::factor(Label, levels = Label))
  
  fill_color <- component_colors[[comp]]
  if (base::is.null(fill_color) || base::is.na(fill_color)) {
    fill_color <- "#999999"
  }
  
  p <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Label, y = Loading)) +
    ggplot2::geom_col(fill = fill_color, alpha = 0.9, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = base::paste0("Top ", top_n, " Contributing Proteins - ", comp),
      x = NULL,
      y = "Loading"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  print(p)
  
  ggplot2::ggsave(
    filename = file.path(output_dir_top, base::paste0("Top20_Contributing_", comp, ".png")),
    plot = p,
    dpi = 300,
    width = 10,
    height = 6,
    bg = "white"
  )
}

cat("\nTop contributing proteins export complete\n")
cat("──────────────────────────────────────────\n")
cat("Components reported: ", base::paste(unique(top_contributing_all$Component), collapse = ", "), "\n", sep = "")
cat("Top proteins per component: ", top_n, "\n", sep = "")
cat("Included original columns from `data`: ", length(all_data_columns), "\n", sep = "")
cat("Combined TSV: ", file.path(output_dir_top, "Top20_Contributing_Proteins_All_Components.tsv"), "\n", sep = "")
cat("Summary TSV: ", file.path(output_dir_top, "Top20_Contributing_Proteins_Summary.tsv"), "\n", sep = "")
cat("Output folder: ", normalizePath(output_dir_top), "\n", sep = "")
cat("Done.\n")

################################################################################

################################################################################
# STEP 3
# Annotation + enrichment for top 20 proteins from sPLS-DA Comp1 only
# P257_09 | HC / BCN / BCP comparison
################################################################################

# ---- [Libraries] ----
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(KEGGREST)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)

# ---- [Working Directory] ----
# Match STEP 1 / STEP 2 location
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Directories] ----
top20_dir <- "Top_Contributing_Proteins_P257_09"
step3_dir <- file.path(top20_dir, "Step3_Annotation_Enrichment_Comp1")
dir.create(step3_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Input Files] ----
top20_combined_file <- file.path(top20_dir, "Top20_Contributing_Proteins_All_Components.tsv")
top20_comp1_file <- file.path(top20_dir, "Top20_Contributing_Comp1.tsv")

# ---- [Helper: annotate species] ----
annotate_species <- function(df) {
  if (!"Protein_IDs" %in% colnames(df)) {
    stop("Column `Protein_IDs` not found in input file.")
  }
  
  df %>%
    dplyr::mutate(
      First_ID = sub(";.*", "", Protein_IDs),
      Species = dplyr::case_when(
        grepl("^hsa:", First_ID, ignore.case = TRUE) ~ "Human",
        grepl("^sav:|^sau:", First_ID, ignore.case = TRUE) ~ "SA",
        
        # Common UniProt-like Human accessions
        grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", First_ID) ~ "Human",
        grepl("^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$", First_ID) ~ "Human",
        grepl("^A0A[0-9A-Z]+$", First_ID) ~ "Human",
        
        # Loose fallback for SA-like IDs
        grepl("^SA", First_ID, ignore.case = TRUE) ~ "SA",
        
        TRUE ~ "Unknown"
      )
    )
}

# ---- [Helper: save annotated table] ----
save_annotated_table <- function(df, outfile) {
  openxlsx::write.xlsx(df, outfile, rowNames = FALSE)
}

# ---- [Helper: enrichment] ----
run_enrichment_from_df <- function(sig_data, out_dir, label_suffix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (nrow(sig_data) == 0) {
    cat("No rows found for:", label_suffix, "\n")
    return(invisible(NULL))
  }
  
  if (!all(c("Protein_IDs", "First_ID", "Species") %in% colnames(sig_data))) {
    stop("Input data must contain Protein_IDs, First_ID, and Species columns.")
  }
  
  # ---- [Split IDs by species] ----
  human_ids_raw <- unique(stats::na.omit(sig_data$First_ID[sig_data$Species == "Human"]))
  sa_ids_raw    <- unique(stats::na.omit(sig_data$First_ID[sig_data$Species == "SA"]))
  
  # Remove prefixes if present
  human_ids <- gsub("^hsa:", "", human_ids_raw, ignore.case = TRUE)
  sa_ids    <- gsub("^sav:|^sau:", "", sa_ids_raw, ignore.case = TRUE)
  
  # ---- [Write input summary] ----
  summary_df <- base::data.frame(
    Category = c("Total_rows", "Human_IDs_raw", "SA_IDs_raw", "Unknown_rows"),
    Count = c(
      nrow(sig_data),
      length(human_ids),
      length(sa_ids),
      sum(sig_data$Species == "Unknown", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  openxlsx::write.xlsx(
    summary_df,
    file = file.path(out_dir, paste0("Input_summary_", label_suffix, ".xlsx")),
    rowNames = FALSE
  )
  
  # --------------------------------------------------------------------------
  # HUMAN ENRICHMENT
  # --------------------------------------------------------------------------
  valid_uniprot <- character(0)
  human_entrez <- NULL
  
  if (length(human_ids) > 0) {
    valid_uniprot <- human_ids[human_ids %in% AnnotationDbi::keys(org.Hs.eg.db, keytype = "UNIPROT")]
    
    if (length(valid_uniprot) > 0) {
      human_entrez <- tryCatch({
        clusterProfiler::bitr(
          valid_uniprot,
          fromType = "UNIPROT",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        )
      }, error = function(e) {
        cat("Human ID conversion failed for", label_suffix, ":", e$message, "\n")
        NULL
      })
      
      if (!is.null(human_entrez) && nrow(human_entrez) > 0) {
        openxlsx::write.xlsx(
          human_entrez,
          file = file.path(out_dir, paste0("Human_ID_conversion_", label_suffix, ".xlsx")),
          rowNames = FALSE
        )
        
        human_genes <- unique(stats::na.omit(human_entrez$ENTREZID))
        
        if (length(human_genes) > 0) {
          # ---- [GO enrichment] ----
          go_enrich <- tryCatch({
            clusterProfiler::enrichGO(
              gene = human_genes,
              OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont = "ALL",
              pAdjustMethod = "BH",
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              readable = TRUE
            )
          }, error = function(e) {
            cat("GO enrichment failed for", label_suffix, ":", e$message, "\n")
            NULL
          })
          
          if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
            openxlsx::write.xlsx(
              as.data.frame(go_enrich),
              file = file.path(out_dir, paste0("GO_Human_enrichment_", label_suffix, ".xlsx")),
              rowNames = FALSE
            )
          }
          
          # ---- [KEGG enrichment: human] ----
          kegg_human <- tryCatch({
            clusterProfiler::enrichKEGG(
              gene = human_genes,
              organism = "hsa",
              keyType = "kegg",
              pvalueCutoff = 1,
              qvalueCutoff = 1
            )
          }, error = function(e) {
            cat("Human KEGG enrichment failed for", label_suffix, ":", e$message, "\n")
            NULL
          })
          
          if (!is.null(kegg_human) && nrow(as.data.frame(kegg_human)) > 0) {
            kegg_human_df <- tryCatch({
              as.data.frame(clusterProfiler::setReadable(
                kegg_human,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID"
              ))
            }, error = function(e) {
              as.data.frame(kegg_human)
            })
            
            openxlsx::write.xlsx(
              kegg_human_df,
              file = file.path(out_dir, paste0("KEGG_Human_enrichment_", label_suffix, ".xlsx")),
              rowNames = FALSE
            )
          }
          
          # ---- [Reactome enrichment: human] ----
          reactome_human <- tryCatch({
            ReactomePA::enrichPathway(
              gene = human_genes,
              organism = "human",
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              readable = TRUE
            )
          }, error = function(e) {
            cat("Reactome enrichment failed for", label_suffix, ":", e$message, "\n")
            NULL
          })
          
          if (!is.null(reactome_human) && nrow(as.data.frame(reactome_human)) > 0) {
            openxlsx::write.xlsx(
              as.data.frame(reactome_human),
              file = file.path(out_dir, paste0("Reactome_Human_enrichment_", label_suffix, ".xlsx")),
              rowNames = FALSE
            )
          }
        }
      } else {
        cat("No human ENTREZ IDs mapped for:", label_suffix, "\n")
      }
    } else {
      cat("No valid human UNIPROT IDs for:", label_suffix, "\n")
    }
  } else {
    cat("No human proteins found for:", label_suffix, "\n")
  }
  
  # --------------------------------------------------------------------------
  # S. AUREUS KEGG ENRICHMENT
  # --------------------------------------------------------------------------
  if (length(sa_ids) > 0) {
    sa_kegg <- tryCatch({
      clusterProfiler::enrichKEGG(
        gene = sa_ids,
        organism = "sau",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
    }, error = function(e) {
      cat("S. aureus KEGG enrichment failed for", label_suffix, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(sa_kegg) && nrow(as.data.frame(sa_kegg)) > 0) {
      openxlsx::write.xlsx(
        as.data.frame(sa_kegg),
        file = file.path(out_dir, paste0("KEGG_SA_enrichment_", label_suffix, ".xlsx")),
        rowNames = FALSE
      )
    }
  } else {
    cat("No S. aureus proteins found for:", label_suffix, "\n")
  }
  
  invisible(NULL)
}

# ---- [Process top 20 proteins from Comp1 only] ----
process_comp1_top20 <- function(top20_combined_file, top20_comp1_file, base_out_dir) {
  
  cat("\n===========================================\n")
  cat("Processing top 20 proteins for sPLS-DA Comp1...\n")
  cat("===========================================\n")
  
  top20_data <- NULL
  
  # Priority 1: specific Comp1 file from STEP 2
  if (file.exists(top20_comp1_file)) {
    top20_data <- utils::read.delim(
      top20_comp1_file,
      header = TRUE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (file.exists(top20_combined_file)) {
    # Fallback: read combined file and filter Comp1
    top20_all <- utils::read.delim(
      top20_combined_file,
      header = TRUE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!"Component" %in% colnames(top20_all)) {
      stop("Column `Component` not found in combined top20 data.")
    }
    
    top20_data <- top20_all %>%
      dplyr::filter(Component == "Comp1")
  }
  
  if (is.null(top20_data) || nrow(top20_data) == 0) {
    cat("No Comp1 top20 input file found or Comp1 table is empty.\n")
    cat("Expected either:\n")
    cat("  ", top20_comp1_file, "\n", sep = "")
    cat("or Comp1 rows inside:\n")
    cat("  ", top20_combined_file, "\n", sep = "")
    return(invisible(NULL))
  }
  
  if (!"Protein_IDs" %in% colnames(top20_data)) {
    stop("Column `Protein_IDs` not found in Comp1 top20 data.")
  }
  
  comp1_annot <- annotate_species(top20_data)
  
  comp1_dir <- file.path(base_out_dir, "Comp1")
  dir.create(comp1_dir, recursive = TRUE, showWarnings = FALSE)
  
  annot_file <- file.path(comp1_dir, "Top20_Comp1_annotated_for_enrichment.xlsx")
  save_annotated_table(comp1_annot, annot_file)
  
  # Save a clean summary table
  summary_cols <- c(
    "Component",
    "Rank",
    "Gene_names",
    "Feature_ID_for_join",
    "Protein_IDs",
    "Majority_protein_IDs",
    "Protein_names",
    "Loading",
    "Abs_Loading",
    "First_ID",
    "Species"
  )
  summary_cols <- summary_cols[summary_cols %in% colnames(comp1_annot)]
  
  comp1_summary <- comp1_annot %>%
    dplyr::select(dplyr::all_of(summary_cols))
  
  openxlsx::write.xlsx(
    comp1_summary,
    file = file.path(comp1_dir, "Top20_Comp1_Annotation_Summary.xlsx"),
    rowNames = FALSE
  )
  
  run_enrichment_from_df(
    sig_data = comp1_annot,
    out_dir = comp1_dir,
    label_suffix = "Top20_Comp1_HC_BCN_BCP"
  )
  
  cat("Processed Comp1 with", nrow(comp1_annot), "proteins\n")
  invisible(NULL)
}

# ---- [Run Comp1 annotation + enrichment] ----
process_comp1_top20(
  top20_combined_file = top20_combined_file,
  top20_comp1_file = top20_comp1_file,
  base_out_dir = step3_dir
)

# ---- [Final report] ----
cat("\n====================================\n")
cat("STEP 3 annotation + enrichment done\n")
cat("====================================\n")
cat("Analysis: Top 20 proteins from sPLS-DA Comp1 only\n")
cat("Comparison: HC / BCN / BCP\n")
cat("Output folder:\n")
cat(normalizePath(step3_dir), "\n")
cat("\nGenerated analyses:\n")
cat("1) Annotation table for Comp1 top 20 proteins\n")
cat("2) Human enrichment: GO, KEGG, Reactome\n")
cat("3) S. aureus enrichment: KEGG\n")
cat("\nAll STEP 3 analyses completed.\n")

################################################################################

################################################################################
# STEP 4
# Fisher's exact test only
# Runs for:
# 1) Top 20 proteins from sPLS-DA Comp1 only
#
# P257_09 | HC / BCN / BCP comparison
#
# NOTE:
# - No shared proteins section
# - odds_ratio removed to avoid Excel #NUM! from infinite estimates
################################################################################

# ---- [Libraries] ----
library(readxl)
library(readr)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Set Paths] ----
base_dir <- getwd()

# STEP 3 Comp1 annotation/enrichment directory
comp1_dir <- file.path(
  base_dir,
  "Top_Contributing_Proteins_P257_09",
  "Step3_Annotation_Enrichment_Comp1",
  "Comp1"
)

# STEP 3 annotated Comp1 file
comp1_sig_file <- file.path(
  comp1_dir,
  "Top20_Comp1_annotated_for_enrichment.xlsx"
)

# STEP 4 output directory
comp1_output_dir <- file.path(
  comp1_dir,
  "Fisher_Test"
)

# Background file from STEP 1
background_file_txt <- file.path(
  base_dir,
  "08.12.25_NEW_Adjusted_Normalized_sepsis_P257_09.txt"
)

dir.create(comp1_output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Load Background Data] ----
if (file.exists(background_file_txt)) {
  background_data <- utils::read.table(
    background_file_txt,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  stop("Background file not found: ", background_file_txt)
}

# ---- [Prepare Background Gene Names] ----
if (!"Gene_names" %in% colnames(background_data)) {
  stop("Column `Gene_names` not found in background data.")
}
if (!"Protein_IDs" %in% colnames(background_data)) {
  stop("Column `Protein_IDs` not found in background data.")
}

background_data$Gene_names <- ifelse(
  is.na(background_data$Gene_names) | trimws(background_data$Gene_names) == "",
  background_data$Protein_IDs,
  background_data$Gene_names
)

background_data$Gene_names <- make.unique(as.character(background_data$Gene_names))

# ---- [Helper: prepare significant genes] ----
prepare_sig_genes <- function(sig_df) {
  
  if (!"Annotated_Gene" %in% colnames(sig_df)) {
    sig_df$Annotated_Gene <- NA_character_
  }
  
  # Prefer Gene_names if present
  if ("Gene_names" %in% colnames(sig_df)) {
    sig_df$Annotated_Gene <- ifelse(
      is.na(sig_df$Annotated_Gene) | trimws(sig_df$Annotated_Gene) == "",
      sig_df$Gene_names,
      sig_df$Annotated_Gene
    )
  }
  
  # Fallback to Protein_IDs if still empty
  if ("Protein_IDs" %in% colnames(sig_df)) {
    sig_df$Annotated_Gene <- ifelse(
      is.na(sig_df$Annotated_Gene) | trimws(sig_df$Annotated_Gene) == "",
      sig_df$Protein_IDs,
      sig_df$Annotated_Gene
    )
  }
  
  sig_df$Annotated_Gene <- trimws(as.character(sig_df$Annotated_Gene))
  sig_df$Annotated_Gene[sig_df$Annotated_Gene == ""] <- NA_character_
  
  unique(na.omit(sig_df$Annotated_Gene))
}

# ---- [Helper: parse enrichment gene list] ----
parse_gene_list <- function(x) {
  vals <- strsplit(as.character(x), "[;/,]")[[1]]
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  unique(vals)
}

# ---- [Core Fisher function] ----
run_fisher_test_set <- function(sig_file,
                                enrichment_dir,
                                output_dir,
                                label_suffix,
                                background_data,
                                sources = c("GO_Human", "KEGG_Human", "Reactome_Human")) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!file.exists(sig_file)) {
    cat("Significant protein file not found:", sig_file, "\n")
    return(invisible(NULL))
  }
  
  sig_df <- readxl::read_xlsx(sig_file)
  sig_genes <- prepare_sig_genes(sig_df)
  
  background_genes <- unique(na.omit(trimws(as.character(background_data$Gene_names))))
  total_background_size <- length(background_genes)
  
  if (length(sig_genes) == 0) {
    cat("No valid significant genes found in:", sig_file, "\n")
    return(invisible(NULL))
  }
  
  if (total_background_size == 0) {
    cat("No valid background genes found.\n")
    return(invisible(NULL))
  }
  
  summary_df <- data.frame(
    Metric = c("Sig_Genes", "Background_Genes"),
    Count = c(length(sig_genes), total_background_size),
    stringsAsFactors = FALSE
  )
  
  openxlsx::write.xlsx(
    summary_df,
    file.path(output_dir, paste0("Fisher_Input_Summary_", label_suffix, ".xlsx")),
    rowNames = FALSE
  )
  
  for (source in sources) {
    
    enrich_file <- file.path(
      enrichment_dir,
      paste0(source, "_enrichment_", label_suffix, ".xlsx")
    )
    
    if (!file.exists(enrich_file)) {
      cat("Skipping missing file:", enrich_file, "\n")
      next
    }
    
    enrich_data <- readxl::read_xlsx(enrich_file)
    
    if (!"geneID" %in% colnames(enrich_data)) {
      cat("No geneID column in", enrich_file, "\n")
      next
    }
    
    fisher_results <- vector("list", nrow(enrich_data))
    keep_idx <- 0
    
    for (i in seq_len(nrow(enrich_data))) {
      
      term <- if ("Description" %in% colnames(enrich_data)) {
        as.character(enrich_data$Description[i])
      } else {
        paste0("Term_", i)
      }
      
      term_genes_raw <- parse_gene_list(enrich_data$geneID[i])
      term_genes <- intersect(term_genes_raw, background_genes)
      sig_in_pathway <- intersect(term_genes, sig_genes)
      
      a <- length(sig_in_pathway)
      b <- length(sig_genes) - a
      c <- length(term_genes) - a
      d <- total_background_size - a - b - c
      
      if (any(c(a, b, c, d) < 0)) next
      
      mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      
      ft <- tryCatch(
        stats::fisher.test(mat),
        error = function(e) NULL
      )
      
      if (is.null(ft)) next
      
      gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
      bg_ratio <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
      enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA_real_)
      
      keep_idx <- keep_idx + 1
      fisher_results[[keep_idx]] <- data.frame(
        Pathway = term,
        Sig_In_Pathway = a,
        Sig_NotIn_Pathway = b,
        Nonsig_In_Pathway = c,
        Nonsig_NotIn_Pathway = d,
        Genes = paste(sig_in_pathway, collapse = ";"),
        p_value = ft$p.value,
        enrichment_factor = enrichment_factor,
        stringsAsFactors = FALSE
      )
    }
    
    fisher_results <- fisher_results[seq_len(keep_idx)]
    
    if (length(fisher_results) > 0) {
      fisher_results_df <- dplyr::bind_rows(fisher_results) %>%
        dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) %>%
        dplyr::arrange(p_adj, p_value)
      
      out_file <- file.path(
        output_dir,
        paste0("Fisher_Pathway_", source, "_", label_suffix, ".xlsx")
      )
      
      openxlsx::write.xlsx(fisher_results_df, out_file, rowNames = FALSE)
      cat("Fisher result saved:", out_file, "\n")
    } else {
      cat("No valid Fisher results for:", source, " / ", label_suffix, "\n")
    }
  }
  
  invisible(NULL)
}

# ---- [Run Fisher test for Comp1 only] ----
run_fisher_test_set(
  sig_file = comp1_sig_file,
  enrichment_dir = comp1_dir,
  output_dir = comp1_output_dir,
  label_suffix = "Top20_Comp1_HC_BCN_BCP",
  background_data = background_data,
  sources = c("GO_Human", "KEGG_Human", "Reactome_Human")
)

# ---- [Final Report] ----
cat("\n====================================\n")
cat("STEP 4 Fisher exact test completed\n")
cat("====================================\n")
cat("Processed:\n")
cat("1) Top 20 proteins from Comp1 only\n")
cat("\nFisher outputs were generated for:\n")
cat("- GO_Human\n")
cat("- KEGG_Human\n")
cat("- Reactome_Human\n")
cat("\nMain output folder:\n")
cat(normalizePath(comp1_output_dir, winslash = "/"), "\n")
cat("\nSTEP 4 done.\n")


################################################################################

################################################################################
# STEP 5
# Enrichment Heatmap Plotting for Fisher's Exact Test Results
# P257_09 | HC / BCN / BCP comparison
# Saves plots inside the Comp1 Fisher_Test folder
# Fixes long-label problem
################################################################################

# ---- [Libraries] ----
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Set Paths] ----
base_dir <- getwd()

comp1_dir <- file.path(
  base_dir,
  "Top_Contributing_Proteins_P257_09",
  "Step3_Annotation_Enrichment_Comp1",
  "Comp1"
)

annot_file <- file.path(
  comp1_dir,
  "Top20_Comp1_annotated_for_enrichment.xlsx"
)

fisher_dir <- file.path(
  comp1_dir,
  "Fisher_Test"
)

# ---- [Settings] ----
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
label_suffix <- "Top20_Comp1_HC_BCN_BCP"

p_adj_cutoff <- 0.05
min_sig_in_pathway <- 3
pathways_per_page <- 20
use_sepsis_filter <- TRUE

sepsis_keywords <- c(
  "sepsis", "infection", "endotoxin", "bacteria",
  "pathogen", "defense", "acute", "positive regulation of innate immune response", 
  "acute inflammation response", "vesicle", "exosome",
  "positive regulation of exocytosis", "extracellular vesicle", 
  "blood coagulation",
  "positive regulation of secrtetion", "blood microparticle", 
  "leukocyte chemotaxis", 
  "lytic vacuole membrane", 
  "secretory granule membrane", "complement activation",
  "positive regulation of adaptive immune response", 
  "defense response to Gram-positive bacterium",
  "positive regulation interleukin-1 production", 
  "macrophage activation", "B cell mediated immunity",
  "toll-like receptor signaling pathway"
)

# ---- [Helper: shorten ugly labels] ----
clean_gene_label <- function(x, max_chars = 30) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  x <- ifelse(
    !is.na(x) & nchar(x) > max_chars,
    paste0(substr(x, 1, max_chars), "..."),
    x
  )
  
  x
}

# ---- [Helper: remove accession-list style labels] ----
is_bad_heatmap_label <- function(x) {
  x <- as.character(x)
  
  too_long <- nchar(x) > 30
  many_ids <- stringr::str_count(x, ";") >= 2
  accession_like <- grepl("^[A-Z0-9.;_-]+$", x)
  
  too_long | (many_ids & accession_like)
}

# ---- [Helper: prepare significant genes from annotated file] ----
get_sig_genes <- function(annot_file) {
  if (!file.exists(annot_file)) {
    cat("Missing annotated file:", annot_file, "\n")
    return(character(0))
  }
  
  annot_df <- readxl::read_xlsx(annot_file)
  
  if (!"Annotated_Gene" %in% colnames(annot_df)) {
    annot_df$Annotated_Gene <- NA_character_
  }
  
  if ("Gene_names" %in% colnames(annot_df)) {
    annot_df$Annotated_Gene <- ifelse(
      is.na(annot_df$Annotated_Gene) | trimws(annot_df$Annotated_Gene) == "",
      annot_df$Gene_names,
      annot_df$Annotated_Gene
    )
  }
  
  if ("Protein_IDs" %in% colnames(annot_df)) {
    annot_df$Annotated_Gene <- ifelse(
      is.na(annot_df$Annotated_Gene) | trimws(annot_df$Annotated_Gene) == "",
      annot_df$Protein_IDs,
      annot_df$Annotated_Gene
    )
  }
  
  sig_genes <- unique(trimws(as.character(annot_df$Annotated_Gene)))
  sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]
  
  # remove labels that break plotting
  sig_genes <- sig_genes[!is_bad_heatmap_label(sig_genes)]
  
  sig_genes
}

# ---- [Helper: detect Fisher count column] ----
detect_sig_count_col <- function(df) {
  candidates <- c("Sig_protein_volcano", "Sig_In_Pathway")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ---- [Helper: load and combine Fisher result files] ----
load_fisher_results <- function(fisher_dir, label_suffix, sources) {
  combined_all <- list()
  
  for (source in sources) {
    fisher_file <- file.path(
      fisher_dir,
      paste0("Fisher_Pathway_", source, "_", label_suffix, ".xlsx")
    )
    
    if (!file.exists(fisher_file)) {
      cat("Missing:", fisher_file, "\n")
      next
    }
    
    fisher_df <- openxlsx::read.xlsx(fisher_file) %>% as.data.frame()
    
    if (nrow(fisher_df) == 0) next
    if (!"Genes" %in% colnames(fisher_df)) next
    if (!"Pathway" %in% colnames(fisher_df)) next
    if (!"p_adj" %in% colnames(fisher_df)) next
    
    sig_count_col <- detect_sig_count_col(fisher_df)
    if (is.na(sig_count_col)) next
    
    fisher_clean <- fisher_df %>%
      dplyr::filter(
        !is.na(p_adj),
        p_adj < p_adj_cutoff,
        .data[[sig_count_col]] >= min_sig_in_pathway,
        !is.na(Genes),
        Genes != ""
      ) %>%
      dplyr::mutate(
        Pathway = stringr::str_trim(Pathway),
        Source = gsub("_Human", "", source),
        Pathway = paste0(stringr::str_wrap(Pathway, 70), " [", Source, "]"),
        log10_p = -log10(p_adj)
      ) %>%
      tidyr::separate_rows(Genes, sep = "[;,/\\s]+") %>%
      dplyr::mutate(Gene = stringr::str_trim(Genes)) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::distinct(Pathway, Gene, log10_p, Source)
    
    combined_all[[source]] <- fisher_clean
  }
  
  dplyr::bind_rows(combined_all)
}

# ---- [Helper: plot paged heatmaps] ----
plot_heatmap_pages <- function(merged_df,
                               sig_genes,
                               dataset_name,
                               save_dir,
                               use_sepsis_filter = TRUE,
                               sepsis_keywords = NULL,
                               pathways_per_page = 20) {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (nrow(merged_df) == 0) {
    cat("No Fisher enrichment rows for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  if (use_sepsis_filter) {
    pathway_keep <- unique(
      merged_df$Pathway[
        stringr::str_detect(
          tolower(merged_df$Pathway),
          paste(tolower(sepsis_keywords), collapse = "|")
        )
      ]
    )
    
    merged_df <- merged_df %>% dplyr::filter(Pathway %in% pathway_keep)
    
    if (nrow(merged_df) == 0) {
      cat("No sepsis-related pathways found for:", dataset_name, "\n")
      return(invisible(NULL))
    }
  }
  
  fisher_genes <- sort(unique(merged_df$Gene))
  sig_genes <- intersect(sig_genes, fisher_genes)
  
  if (length(sig_genes) == 0) {
    cat("No overlapping genes between annotation and Fisher results for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  gene_label_map <- data.frame(
    Gene = sig_genes,
    Gene_Label = clean_gene_label(sig_genes, max_chars = 25),
    stringsAsFactors = FALSE
  )
  
  gene_label_map$Gene_Label <- make.unique(gene_label_map$Gene_Label)
  
  pathway_levels <- unique(merged_df$Pathway)
  
  full_grid <- expand.grid(
    Gene = sig_genes,
    Pathway = pathway_levels,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(merged_df, by = c("Gene", "Pathway")) %>%
    dplyr::left_join(gene_label_map, by = "Gene")
  
  if (nrow(full_grid) == 0) {
    cat("Nothing to plot for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  pathway_ranking <- full_grid %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarize(avg_log10_p = mean(log10_p, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(avg_log10_p))
  
  if (nrow(pathway_ranking) == 0) {
    cat("No ranked pathways for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  gene_levels <- gene_label_map$Gene_Label[match(sort(sig_genes), gene_label_map$Gene)]
  
  full_grid <- full_grid %>%
    dplyr::mutate(
      Pathway = factor(Pathway, levels = pathway_ranking$Pathway),
      Gene_Label = factor(Gene_Label, levels = gene_levels)
    )
  
  pathway_chunks <- split(
    pathway_ranking$Pathway,
    ceiling(seq_along(pathway_ranking$Pathway) / pathways_per_page)
  )
  
  for (i in seq_along(pathway_chunks)) {
    current_paths <- pathway_chunks[[i]]
    
    chunk_data <- full_grid %>%
      dplyr::filter(Pathway %in% current_paths) %>%
      dplyr::mutate(Pathway = factor(Pathway, levels = rev(current_paths)))
    
    plot_file <- file.path(
      save_dir,
      paste0(
        gsub("[^A-Za-z0-9_]+", "_", dataset_name),
        "_Heatmap_Page",
        i,
        ".png"
      )
    )
    
    p <- ggplot2::ggplot(chunk_data, ggplot2::aes(x = Gene_Label, y = Pathway, fill = log10_p)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
      ggplot2::scale_fill_gradient(
        low = "lightblue",
        high = "red",
        na.value = "black",
        name = "-log10(p.adj)"
      ) +
      ggplot2::labs(
        title = NULL,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 5) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 32),
        axis.text.y = ggplot2::element_text(size = 32),
        axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = ggplot2::element_text(size = 32, margin = ggplot2::margin(b = 15)),
        legend.text = ggplot2::element_text(size = 32),
        legend.key.height = grid::unit(2, "cm"),
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )
    
    ggplot2::ggsave(
      filename = plot_file,
      plot = p,
      width = max(30, length(unique(chunk_data$Gene_Label)) * 1),
      height = max(15, length(current_paths) * 0.8),
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    
    cat("Saved:", plot_file, "\n")
  }
  
  invisible(NULL)
}

# ---- [Run STEP 5 for Comp1 only] ----
comp1_sig_genes <- get_sig_genes(annot_file)

if (length(comp1_sig_genes) == 0) {
  cat("No significant genes found for Comp1.\n")
} else if (!dir.exists(fisher_dir)) {
  cat("Fisher_Test folder not found:", fisher_dir, "\n")
} else {
  
  comp1_merged <- load_fisher_results(
    fisher_dir = fisher_dir,
    label_suffix = label_suffix,
    sources = sources
  )
  
  plot_heatmap_pages(
    merged_df = comp1_merged,
    sig_genes = comp1_sig_genes,
    dataset_name = "Top20_Comp1_HC_BCN_BCP",
    save_dir = fisher_dir,
    use_sepsis_filter = use_sepsis_filter,
    sepsis_keywords = sepsis_keywords,
    pathways_per_page = pathways_per_page
  )
}

# ---- [Final Report] ----
cat("\n=================================================\n")
cat("STEP 5 enrichment heatmap plotting completed\n")
cat("=================================================\n")
cat("Processed dataset: Top 20 proteins from Comp1 only\n")
cat("Comparison: HC / BCN / BCP\n")
cat("Plots saved inside:\n")
cat(normalizePath(fisher_dir, winslash = "/"), "\n")
cat("Long accession-style labels were filtered/shortened for plotting.\n")
cat("\nSTEP 5 done.\n")
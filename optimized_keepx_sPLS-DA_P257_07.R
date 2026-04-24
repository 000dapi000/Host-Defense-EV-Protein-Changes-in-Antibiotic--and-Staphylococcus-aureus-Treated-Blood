library(mixOmics)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(BiocParallel)

# ---- [Global Reproducibility Seed] ----
SEED <- 1
set.seed(SEED)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis")

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Load Expression Data] ----
data <- utils::read.table(
  file = "Adjusted_Normalized_P257_07.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- [Check Required Columns] ----
required_cols <- c("Gene_names", "Protein_IDs")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# ---- [Fix Gene Names] ----
data$Gene_names <- ifelse(
  is.na(data$Gene_names) | trimws(data$Gene_names) == "",
  data$Protein_IDs,
  data$Gene_names
)

data$Gene_names <- make.unique(as.character(data$Gene_names))
rownames(data) <- data$Gene_names

# ---- [Define Groups] ----
group_cols <- list(
  MOCK         = grep("^MOCK_[1-6]$", colnames(data), value = TRUE),
  SA           = grep("^SA_[1-6]$", colnames(data), value = TRUE),
  PipTazo      = grep("^PipTazo_[1-6]$", colnames(data), value = TRUE),
  Vancomycin   = grep("^Vancomycin_[1-6]$", colnames(data), value = TRUE),
  Moxifloxacin = grep("^Moxifloxacin_[1-6]$", colnames(data), value = TRUE)
)

selected_groups <- names(group_cols)

# ---- [Check Group Sizes] ----
group_sizes <- vapply(group_cols, length, integer(1))
print(group_sizes)

if (any(group_sizes == 0)) {
  stop("One or more groups have zero matched columns. Please check sample names.")
}

# ---- [Create X and Y] ----
sample_columns <- unlist(group_cols, use.names = FALSE)
X_data <- data[, sample_columns, drop = FALSE]
X <- as.data.frame(t(X_data), stringsAsFactors = FALSE)
sample_names <- rownames(X)

Y <- dplyr::case_when(
  grepl("^MOCK", sample_names) ~ "MOCK",
  grepl("^SA", sample_names) ~ "SA",
  grepl("^PipTazo", sample_names) ~ "PipTazo",
  grepl("^Vancomycin", sample_names) ~ "Vancomycin",
  grepl("^Moxifloxacin", sample_names) ~ "Moxifloxacin",
  TRUE ~ NA_character_
)

valid_idx <- which(!is.na(Y))
if (length(valid_idx) == 0) {
  stop("No valid samples matched the defined groups.")
}

X <- X[valid_idx, , drop = FALSE]
sample_names <- sample_names[valid_idx]
Y <- factor(Y[valid_idx], levels = selected_groups)

# ---- [Ensure Numeric Data for mixOmics] ----
X <- as.data.frame(lapply(X, as.numeric), stringsAsFactors = FALSE)
rownames(X) <- sample_names

# ---- [Remove Features with NA or Zero Variance] ----
valid_features <- vapply(
  X,
  function(z) {
    all(is.finite(z)) && stats::sd(z, na.rm = TRUE) > 0
  },
  logical(1)
)

X <- X[, valid_features, drop = FALSE]

if (ncol(X) == 0) {
  stop("No valid variables remain after filtering for finite values and non-zero variance.")
}

cat("Samples:", nrow(X), "\n")
cat("Variables after filtering:", ncol(X), "\n")
cat("Class counts:\n")
print(table(Y))

# ---- [Color Palette] ----
group_colors <- c(
  "MOCK"         = "#1b9e77",
  "SA"           = "#d95f02",
  "PipTazo"      = "#7570b3",
  "Vancomycin"   = "#e7298a",
  "Moxifloxacin" = "#66a61e"
)

# ---- [Custom Theme] ----
custom_theme_splsda <- function(base_size = 25) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.position = "none"
    )
}

# ---- [Tuning Parameters] ----
optimal_ncomp <- 4
list_keepX <- c(25, 50, 100)

# ---- [Reproducible Serial Backend] ----
bp_param <- SerialParam()

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

tune_formals <- names(formals(mixOmics::tune.splsda))

if ("seed" %in% tune_formals) {
  tune_args$seed <- SEED
}

if ("BPPARAM" %in% tune_formals) {
  tune_args$BPPARAM <- bp_param
}

tune_result <- do.call(mixOmics::tune.splsda, tune_args)

optimal_keepX <- tune_result$choice.keepX[seq_len(optimal_ncomp)]
print(optimal_keepX)

# ---- [Save Tuning Summary] ----
tuning_summary <- data.frame(
  Component = paste0("Comp", seq_along(optimal_keepX)),
  keepX = as.numeric(optimal_keepX),
  stringsAsFactors = FALSE
)

utils::write.table(
  x = tuning_summary,
  file = file.path(output_dir, "optimal_keepX_summary_P257_07.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---- [keepX Barplot] ----
keepx_df <- data.frame(
  Component = factor(
    paste0("Comp", seq_along(optimal_keepX)),
    levels = paste0("Comp", seq_along(optimal_keepX))
  ),
  keepX = as.numeric(optimal_keepX),
  stringsAsFactors = FALSE
)

gg_keepx <- ggplot(
  keepx_df,
  aes(x = Component, y = keepX)
) +
  geom_col(fill = "#377EB8", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 5) +
  labs(
    title = "Optimal keepX per Component",
    x = "Component",
    y = "Variables Selected"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA)
  )

print(gg_keepx)

ggsave(
  filename = file.path(output_dir, "keepX_per_component_P257_07.png"),
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

perf_formals <- names(formals(mixOmics::perf))

if ("seed" %in% perf_formals) {
  perf_args$seed <- SEED
}

if ("BPPARAM" %in% perf_formals) {
  perf_args$BPPARAM <- bp_param
}

perf_result <- do.call(mixOmics::perf, perf_args)

# ---- [NEW: Extract BER summary] ----
# mixOmics stores error rates by distance; for classification this usually includes
# $error.rate$BER or $error.rate$centroids.dist depending on version/output structure.

ber_by_comp <- NULL

if (!is.null(perf_result$error.rate$BER)) {
  ber_by_comp <- perf_result$error.rate$BER
} else if (!is.null(perf_result$error.rate$centroids.dist)) {
  ber_by_comp <- perf_result$error.rate$centroids.dist
}

if (is.null(ber_by_comp)) {
  warning("Could not extract BER directly from perf_result$error.rate.")
} else {
  ber_summary <- data.frame(
    Component = paste0("Comp", seq_along(ber_by_comp)),
    BER = as.numeric(ber_by_comp),
    stringsAsFactors = FALSE
  )
  
  utils::write.table(
    x = ber_summary,
    file = file.path(output_dir, "BER_by_component_P257_07.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  print(ber_summary)
}

# ---- [NEW: Confusion Matrix from Final Number of Components] ----
conf_mat <- NULL

if (!is.null(perf_result$confusion)) {
  # Some mixOmics versions store confusion matrices indexed by distance and component
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
  conf_mat_df <- as.data.frame.matrix(conf_mat)
  
  utils::write.table(
    x = cbind(True_Class = rownames(conf_mat_df), conf_mat_df),
    file = file.path(output_dir, "confusion_matrix_P257_07.tsv"),
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
  gg_ber <- ggplot(ber_summary, aes(x = Component, y = BER)) +
    geom_col(fill = "#D95F02", alpha = 0.85, width = 0.6) +
    geom_text(aes(label = round(BER, 3)), vjust = -0.5, size = 5) +
    labs(
      title = "Cross-validated BER per Component",
      x = "Component",
      y = "Balanced Error Rate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA)
    )
  
  print(gg_ber)
  
  ggsave(
    filename = file.path(output_dir, "BER_per_component_P257_07.png"),
    plot = gg_ber,
    dpi = 300,
    width = 6,
    height = 5,
    bg = "white"
  )
}

# ---- [Explained Variance Proxy] ----
expl_var <- apply(splsda_model$variates$X^2, 2, sum) / sum(splsda_model$variates$X^2)

# ---- [Robust Ellipse Function] ----
desired_ellipse_level <- 0.95

compute_ellipse <- function(mean_vec, cov_mat, level = 0.95, npoints = 200) {
  cov_mat <- as.matrix(cov_mat)
  
  if (any(!is.finite(cov_mat))) {
    return(NULL)
  }
  
  eps <- 1e-8
  cov_mat <- cov_mat + diag(eps, nrow(cov_mat))
  
  eig <- eigen(cov_mat)
  eig$values[eig$values < eps] <- eps
  
  angles <- seq(0, 2 * pi, length.out = npoints)
  circle <- rbind(cos(angles), sin(angles))
  radius <- sqrt(stats::qchisq(level, df = 2))
  
  transform_mat <- eig$vectors %*% diag(sqrt(eig$values))
  ellipse <- t(radius * transform_mat %*% circle)
  ellipse <- sweep(ellipse, 2, mean_vec, "+")
  
  ellipse_df <- as.data.frame(ellipse, stringsAsFactors = FALSE)
  colnames(ellipse_df) <- c("comp1", "comp2")
  ellipse_df
}

# ---- [Plot All Component Pairs] ----
for (i in seq_len(optimal_ncomp - 1L)) {
  for (j in seq.int(i + 1L, optimal_ncomp)) {
    
    plot_data <- data.frame(
      comp1 = splsda_model$variates$X[, i],
      comp2 = splsda_model$variates$X[, j],
      Group = Y,
      Sample = rownames(X),
      stringsAsFactors = FALSE
    )
    
    group_centroids <- plot_data %>%
      group_by(Group) %>%
      summarise(
        comp1 = mean(comp1, na.rm = TRUE),
        comp2 = mean(comp2, na.rm = TRUE),
        count = n(),
        .groups = "drop"
      )
    
    ellipse_list <- lapply(split(plot_data, plot_data$Group), function(df_group) {
      if (nrow(df_group) < 3) {
        return(NULL)
      }
      
      group_data <- df_group[, c("comp1", "comp2"), drop = FALSE]
      cov_mat <- stats::cov(group_data)
      
      ell <- tryCatch(
        compute_ellipse(
          mean_vec = colMeans(group_data),
          cov_mat = cov_mat,
          level = desired_ellipse_level
        ),
        error = function(e) NULL
      )
      
      if (is.null(ell)) {
        return(NULL)
      }
      
      ell$Group <- unique(df_group$Group)
      ell
    })
    
    ellipse_list <- Filter(Negate(is.null), ellipse_list)
    
    ellipse_data <- if (length(ellipse_list) > 0) {
      bind_rows(ellipse_list)
    } else {
      data.frame(
        comp1 = numeric(0),
        comp2 = numeric(0),
        Group = character(0),
        stringsAsFactors = FALSE
      )
    }
    
    x_lab <- paste0(
      "Component ", i, " (",
      round(expl_var[i] * 100, 1), "%)"
    )
    y_lab <- paste0(
      "Component ", j, " (",
      round(expl_var[j] * 100, 1), "%)"
    )
    
    p <- ggplot(
      plot_data,
      aes(x = comp1, y = comp2, color = Group, fill = Group)
    ) +
      geom_point(size = 4, alpha = 0.9)
    
    if (nrow(ellipse_data) > 0) {
      p <- p +
        geom_polygon(
          data = ellipse_data,
          mapping = aes(
            x = comp1,
            y = comp2,
            group = Group,
            fill = Group
          ),
          alpha = 0.2,
          color = NA,
          inherit.aes = FALSE
        ) +
        geom_path(
          data = ellipse_data,
          mapping = aes(
            x = comp1,
            y = comp2,
            color = Group,
            group = Group
          ),
          linewidth = 1,
          inherit.aes = FALSE
        )
    }
    
    p <- p +
      ggrepel::geom_text_repel(
        data = group_centroids,
        mapping = aes(
          x = comp1,
          y = comp2,
          label = paste0(Group, "\n(n=", count, ")")
        ),
        color = "black",
        size = 10,
        fontface = "bold",
        max.overlaps = 100,
        box.padding = 0.6,
        point.padding = 0.6,
        seed = SEED,
        inherit.aes = FALSE
      ) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      labs(
        title = paste0("sPLS-DA: P257_07 (Comp ", i, " vs Comp ", j, ")"),
        x = x_lab,
        y = y_lab
      ) +
      custom_theme_splsda()
    
    print(p)
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0("sPLS-DA_P257_07_Comp", i, "_vs_Comp", j, ".png")
      ),
      plot = p,
      dpi = 300,
      width = 10,
      height = 10,
      bg = "white"
    )
  }
}

# ---- [Session/Version Log for Reproducibility] ----
sink(file.path(output_dir, "sessionInfo_P257_07.txt"))
cat("Seed used:", SEED, "\n\n")
print(sessionInfo())
sink()

# ---- [Summary Report] ----
cat("\nsPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("Seed Used:", SEED, "\n")
cat("Number of Components Used:", optimal_ncomp, "\n")
cat("keepX per Component:", paste(optimal_keepX, collapse = ", "), "\n")
cat("Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")

if (exists("ber_summary")) {
  cat("Cross-validated BER by component:\n")
  print(ber_summary)
}

n_to_report <- min(2, length(expl_var))
cat(
  "Explained Variance (First ", n_to_report, " components): ",
  paste0(round(expl_var[seq_len(n_to_report)] * 100, 1), collapse = "%, "),
  "%\n",
  sep = ""
)

cat("Output folder:", normalizePath(output_dir), "\n")
cat("\nAll sPLS-DA 2D plots and validation summaries saved successfully.\n")




################################################################################

################################################################################
# STEP 2
# Export top contributing proteins per sPLS-DA component to TSV
# Correct logic:
#   1) use Feature_ID_for_join only to map loadings back to original `data`
#   2) then use Protein_IDs from matched rows to recover annotation gene names
################################################################################

# ---- [Optional Package Checks] ----
required_pkgs_step2 <- c("dplyr", "ggplot2")
missing_pkgs_step2 <- required_pkgs_step2[
  !vapply(required_pkgs_step2, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs_step2) > 0) {
  stop(
    paste0(
      "Please install the following packages before running STEP 2: ",
      paste(missing_pkgs_step2, collapse = ", ")
    )
  )
}

# ---- [Set Parameters] ----
top_n <- 20
n_components_to_report <- min(4, ncol(splsda_model$loadings$X))
output_dir <- "Top_Contributing_Proteins"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotation_file <- "Annotations(Human+all S. aureus).txt"

# ---- [Check Required Objects] ----
if (!exists("data")) {
  stop("Object `data` not found. Please run STEP 1 first.")
}
if (!exists("splsda_model")) {
  stop("Object `splsda_model` not found. Please run STEP 1 first.")
}
if (!"Protein_IDs" %in% colnames(data)) {
  stop("Column `Protein_IDs` not found in `data`.")
}

# ---- [Load Annotation File] ----
Annotation_Human_SA <- utils::read.table(
  annotation_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))

# adjust this if your exact cleaned column name differs
entry_col <- "Entry"
gene_col_candidates <- c("Gene_Names__primary_", "Gene_Names_primary_", "Gene_Names")

if (!entry_col %in% colnames(Annotation_Human_SA)) {
  stop("Column `Entry` not found in annotation file.")
}

gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(Annotation_Human_SA)][1]

if (is.na(gene_col)) {
  stop("Could not find a gene name column in annotation file.")
}

annotation_map <- Annotation_Human_SA |>
  dplyr::transmute(
    Entry = trimws(as.character(.data[[entry_col]])),
    Annot_Gene = trimws(as.character(.data[[gene_col]]))
  ) |>
  dplyr::filter(
    !is.na(.data$Entry),
    .data$Entry != ""
  ) |>
  dplyr::distinct(.data$Entry, .keep_all = TRUE)

# ---- [Helper: resolve Protein_IDs using annotation Entry] ----
resolve_annotation_label <- function(protein_ids_string, annotation_df) {
  if (is.na(protein_ids_string) || trimws(protein_ids_string) == "") {
    return(NA_character_)
  }
  
  ids <- unlist(strsplit(as.character(protein_ids_string), ";", fixed = TRUE))
  ids <- trimws(ids)
  ids <- ids[ids != ""]
  ids <- unique(ids)
  
  if (length(ids) == 0) {
    return(NA_character_)
  }
  
  matched <- annotation_df$Annot_Gene[match(ids, annotation_df$Entry)]
  matched <- matched[!is.na(matched)]
  matched <- trimws(matched)
  matched <- matched[matched != ""]
  matched <- unique(matched)
  
  if (length(matched) == 0) {
    return(NA_character_)
  }
  
  paste(matched, collapse = ";")
}

# ---- [Prepare Original Data for Join Back from Model] ----
data_export <- data

# keep as character
data_export$Protein_IDs <- as.character(data_export$Protein_IDs)

if (!"Gene_names" %in% colnames(data_export)) {
  data_export$Gene_names <- NA_character_
} else {
  data_export$Gene_names <- as.character(data_export$Gene_names)
}

# your original logic from STEP 1
data_export$Gene_names <- ifelse(
  is.na(data_export$Gene_names) | trimws(data_export$Gene_names) == "",
  data_export$Protein_IDs,
  data_export$Gene_names
)

data_export$Gene_names <- make.unique(as.character(data_export$Gene_names))

# IMPORTANT:
# this is ONLY for matching model feature names back to rows in data
data_export$Feature_ID_for_join <- make.names(data_export$Gene_names, unique = TRUE)

# ---- [After row recovery, resolve annotation from Protein_IDs] ----
data_export$Annot_Gene_from_ProteinIDs <- vapply(
  data_export$Protein_IDs,
  resolve_annotation_label,
  character(1),
  annotation_df = annotation_map
)

# choose the final display label
# if Gene_names is just the raw accession string, prefer resolved annotation label
data_export$Gene_names_for_plot <- dplyr::case_when(
  !is.na(data_export$Gene_names) &
    trimws(data_export$Gene_names) != "" &
    data_export$Gene_names != data_export$Protein_IDs ~ data_export$Gene_names,
  
  !is.na(data_export$Annot_Gene_from_ProteinIDs) &
    trimws(data_export$Annot_Gene_from_ProteinIDs) != "" ~ data_export$Annot_Gene_from_ProteinIDs,
  
  TRUE ~ data_export$Protein_IDs
)

# original column order from input
all_data_columns <- colnames(data)

# ---- [Extract Top Features Per Component] ----
top_tables <- lapply(seq_len(n_components_to_report), function(comp) {
  
  loadings_vec <- splsda_model$loadings$X[, comp]
  loadings_vec <- loadings_vec[!is.na(loadings_vec)]
  
  ord <- order(abs(loadings_vec), decreasing = TRUE)
  top_features <- loadings_vec[ord][seq_len(min(top_n, length(loadings_vec)))]
  
  feature_df <- data.frame(
    Feature_ID_for_join = names(top_features),
    Loading = as.numeric(top_features),
    Abs_Loading = abs(as.numeric(top_features)),
    Rank = seq_along(top_features),
    Component = paste0("Comp", comp),
    stringsAsFactors = FALSE
  )
  
  # THIS join is model feature -> original data row
  merged_df <- dplyr::left_join(feature_df, data_export, by = "Feature_ID_for_join")
  
  merged_df <- dplyr::select(
    merged_df,
    .data$Component,
    .data$Rank,
    .data$Gene_names_for_plot,
    .data$Gene_names,
    .data$Annot_Gene_from_ProteinIDs,
    .data$Feature_ID_for_join,
    .data$Loading,
    .data$Abs_Loading,
    dplyr::all_of(all_data_columns),
    dplyr::everything()
  )
  
  merged_df
})

# ---- [Combine All Components] ----
top_contributing_all <- dplyr::bind_rows(top_tables)

# ---- [Check unmatched rows] ----
protein_ids_exists <- "Protein_IDs" %in% colnames(top_contributing_all)
protein_names_exists <- "Protein_names" %in% colnames(top_contributing_all)

if (protein_ids_exists && protein_names_exists) {
  unmatched_df <- dplyr::filter(
    top_contributing_all,
    is.na(.data$Protein_IDs) & is.na(.data$Protein_names)
  )
} else if (protein_ids_exists) {
  unmatched_df <- dplyr::filter(
    top_contributing_all,
    is.na(.data$Protein_IDs)
  )
} else {
  unmatched_df <- top_contributing_all[0, , drop = FALSE]
}

if (nrow(unmatched_df) > 0) {
  utils::write.table(
    x = unmatched_df,
    file = file.path(output_dir, "Unmatched_Top_Features.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  cat("Warning:", nrow(unmatched_df), "top features could not be matched back to `data`.\n")
}

# ---- [Write combined TSV] ----
utils::write.table(
  x = top_contributing_all,
  file = file.path(output_dir, "Top20_Contributing_Proteins_All_Components.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# ---- [Write per-component TSV] ----
for (comp in unique(top_contributing_all$Component)) {
  comp_df <- dplyr::filter(top_contributing_all, .data$Component == comp)
  
  utils::write.table(
    x = comp_df,
    file = file.path(output_dir, paste0("Top20_Contributing_", comp, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

# ---- [Write summary TSV] ----
summary_cols <- c(
  "Component",
  "Rank",
  "Gene_names_for_plot",
  "Gene_names",
  "Annot_Gene_from_ProteinIDs",
  "Feature_ID_for_join",
  "Protein_IDs",
  "Majority_protein_IDs",
  "Protein_names",
  "Loading",
  "Abs_Loading"
)

summary_cols <- summary_cols[summary_cols %in% colnames(top_contributing_all)]

summary_df <- dplyr::select(
  top_contributing_all,
  dplyr::all_of(summary_cols)
)

utils::write.table(
  x = summary_df,
  file = file.path(output_dir, "Top20_Contributing_Proteins_Summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# ---- [Plots] ----
component_colors <- c(
  "Comp1" = "#F8766D",
  "Comp2" = "#7CAE00",
  "Comp3" = "#00BFC4",
  "Comp4" = "#C77CFF",
  "Comp5" = "#FF61C3",
  "Comp6" = "#00BA38"
)

plot_df <- dplyr::mutate(
  top_contributing_all,
  Label = dplyr::case_when(
    !is.na(.data$Gene_names_for_plot) &
      trimws(.data$Gene_names_for_plot) != "" ~ .data$Gene_names_for_plot,
    !is.na(.data$Annot_Gene_from_ProteinIDs) &
      trimws(.data$Annot_Gene_from_ProteinIDs) != "" ~ .data$Annot_Gene_from_ProteinIDs,
    !is.na(.data$Gene_names) &
      trimws(.data$Gene_names) != "" ~ .data$Gene_names,
    !is.na(.data$Protein_IDs) &
      trimws(.data$Protein_IDs) != "" ~ .data$Protein_IDs,
    TRUE ~ .data$Feature_ID_for_join
  )
)

for (comp in unique(plot_df$Component)) {
  comp_df <- plot_df |>
    dplyr::filter(.data$Component == comp) |>
    dplyr::arrange(.data$Abs_Loading) |>
    dplyr::mutate(
      Label_short = ifelse(
        nchar(.data$Label) > 60,
        paste0(substr(.data$Label, 1, 57), "..."),
        .data$Label
      ),
      Label_short = factor(.data$Label_short, levels = .data$Label_short)
    )
  
  fill_col <- if (comp %in% names(component_colors)) component_colors[[comp]] else "#4DAF4A"
  
  p <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Label_short, y = Loading)) +
    ggplot2::geom_col(fill = fill_col, alpha = 0.9, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", top_n, " Contributing Proteins - ", comp),
      x = NULL,
      y = "Loading value"
    ) +
    ggplot2::theme_minimal(base_size = 25) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 25),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 25),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  print(p)
  
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0("Top20_Contributing_", comp, ".png")),
    plot = p,
    dpi = 300,
    width = 12,
    height = 8,
    bg = "white"
  )
}

# ---- [Console Report] ----
cat("\nTop contributing proteins export complete\n")
cat("──────────────────────────────────────────\n")
cat("Components reported: ", paste(unique(top_contributing_all$Component), collapse = ", "), "\n", sep = "")
cat("Top proteins per component: ", top_n, "\n", sep = "")
cat("Included original columns from `data`: ", length(all_data_columns), "\n", sep = "")
cat("Combined TSV: ", file.path(output_dir, "Top20_Contributing_Proteins_All_Components.tsv"), "\n", sep = "")
cat("Summary TSV: ", file.path(output_dir, "Top20_Contributing_Proteins_Summary.tsv"), "\n", sep = "")
cat("Output folder: ", normalizePath(output_dir), "\n", sep = "")
cat("Done.\n")

################################################################################

################################################################################
# STEP 3
# Annotate top Comp1 contributors, annotate SA vs MOCK up-regulated proteins,
# export shared proteins with ALL columns from `data`,
# and create a corrected Venn diagram of overlap
################################################################################

# ---- [Load Libraries] ----
library(mixOmics)
library(dplyr)
library(readr)
library(readxl)
library(openxlsx)
library(VennDiagram)
library(grid)

# ---- [Working Directory & Output Folder] ----
setwd("C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis")
output_dir <- "sPLSDA_2D_Plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Check Required Objects] ----
if (!exists("splsda_model")) {
  stop("Object `splsda_model` not found. Please run STEP 1 first.")
}

if (!exists("data")) {
  stop("Object `data` not found. Please run STEP 1 first.")
}

if (is.null(splsda_model$loadings$X) || ncol(splsda_model$loadings$X) < 1) {
  stop("`splsda_model$loadings$X` is missing or has no components.")
}

if (!"Protein_IDs" %in% colnames(data)) {
  stop("Column `Protein_IDs` not found in `data`.")
}

# ---- [Helper: Clean Column Names Robustly] ----
clean_names <- function(x) {
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[[:punct:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---- [Helper: split semicolon-separated values] ----
split_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(x) == "") {
    return(character(0))
  }
  ids <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  ids <- trimws(ids)
  ids <- ids[ids != ""]
  unique(ids)
}

# ---- [Helper: collapse unique non-empty strings] ----
collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

# ---- [Helper: normalize labels for overlap/Venn] ----
normalize_labels <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  x <- gsub("[[:space:]]+", "", x)
  unique(x)
}

# ---- [Helper: choose best label from a row] ----
pick_best_label <- function(gene_names_for_plot = NA_character_,
                            annotated_gene = NA_character_,
                            gene_names = NA_character_,
                            protein_ids = NA_character_,
                            first_id = NA_character_) {
  vals <- c(gene_names_for_plot, annotated_gene, gene_names, protein_ids, first_id)
  vals <- as.character(vals)
  vals <- trimws(vals)
  vals <- vals[!is.na(vals) & vals != ""]
  if (length(vals) == 0) return(NA_character_)
  vals[1]
}

# ---- [Load & Sanitize Global Annotation Table] ----
Annotation_Human_SA <- utils::read.table(
  file = "Annotations(Human+all S. aureus).txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(Annotation_Human_SA) <- clean_names(colnames(Annotation_Human_SA))

# ---- [Check Required Annotation Columns] ----
required_annotation_cols <- clean_names(c(
  "Entry",
  "Gene Names (primary)",
  "Gene Ontology (GO)",
  "KEGG",
  "Reactome"
))

missing_annotation_cols <- required_annotation_cols[
  !required_annotation_cols %in% colnames(Annotation_Human_SA)
]

if (length(missing_annotation_cols) > 0) {
  cat("Available columns in annotation file:\n")
  print(colnames(Annotation_Human_SA))
  stop(
    paste0(
      "Missing required annotation columns: ",
      paste(missing_annotation_cols, collapse = ", ")
    )
  )
}

# ---- [Build a Clean Annotation Lookup] ----
annot_lookup <- Annotation_Human_SA |>
  dplyr::select(
    .data$Entry,
    .data$Gene_Names_primary,
    .data$Gene_Ontology_GO,
    .data$KEGG,
    .data$Reactome
  ) |>
  dplyr::rename(
    First_ID = .data$Entry,
    Annotated_Gene = .data$Gene_Names_primary,
    Human_GO = .data$Gene_Ontology_GO,
    Human_KEGG = .data$KEGG,
    Human_Reactome = .data$Reactome
  ) |>
  dplyr::mutate(
    First_ID = trimws(as.character(.data$First_ID)),
    Annotated_Gene = trimws(as.character(.data$Annotated_Gene))
  ) |>
  dplyr::distinct(.data$First_ID, .keep_all = TRUE)

# ---- [Prepare original data with ALL columns preserved] ----
all_data_columns <- colnames(data)

data_full <- data
data_full$Protein_IDs <- as.character(data_full$Protein_IDs)

if ("Gene_names" %in% colnames(data_full)) {
  data_full$Gene_names <- as.character(data_full$Gene_names)
} else {
  data_full$Gene_names <- NA_character_
}

data_full$Gene_names <- ifelse(
  is.na(data_full$Gene_names) | trimws(data_full$Gene_names) == "",
  data_full$Protein_IDs,
  data_full$Gene_names
)

# ---- [Expand original data so each semicolon-separated accession gets its own join key] ----
data_expanded <- do.call(
  rbind,
  lapply(seq_len(nrow(data_full)), function(i) {
    ids <- split_ids(data_full$Protein_IDs[i])
    if (length(ids) == 0) {
      ids <- NA_character_
    }
    out <- data_full[rep(i, length(ids)), , drop = FALSE]
    out$First_ID <- ids
    out
  })
)

rownames(data_expanded) <- NULL

# ---- [1) Export Top-20 Comp1 contributors from loadings] ----
df_ld <- data.frame(
  raw_ID = names(splsda_model$loadings$X[, 1]),
  loading = abs(splsda_model$loadings$X[, 1]),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    First_ID = sub(";.*", "", .data$raw_ID)
  )

top20_df <- df_ld |>
  dplyr::group_by(.data$First_ID) |>
  dplyr::slice_max(order_by = .data$loading, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(dplyr::desc(.data$loading)) |>
  dplyr::slice_head(n = 20)

top20_annot <- top20_df |>
  dplyr::left_join(annot_lookup, by = "First_ID") |>
  dplyr::mutate(
    Annotated_Gene = ifelse(
      is.na(.data$Annotated_Gene) | .data$Annotated_Gene == "",
      .data$First_ID,
      .data$Annotated_Gene
    )
  )

utils::write.table(
  x = top20_df,
  file = file.path(output_dir, "Top20_Comp1.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

openxlsx::write.xlsx(
  x = top20_annot,
  file = file.path(output_dir, "Top20_Comp1_annotated.xlsx"),
  rowNames = FALSE
)

# ---- [2) Annotate SA vs MOCK Up-regulated List] ----
sa_up_file <- file.path("SA_vs_MOCK", "Upregulated_in_SA_SA_vs_MOCK.txt")

if (!file.exists(sa_up_file)) {
  stop(
    paste0(
      "File not found: ",
      normalizePath(dirname(sa_up_file), winslash = "/", mustWork = FALSE),
      "/",
      basename(sa_up_file)
    )
  )
}

sa_up_raw <- readr::read_tsv(
  file = sa_up_file,
  show_col_types = FALSE
)

if (!"Protein_IDs" %in% colnames(sa_up_raw)) {
  stop("Column `Protein_IDs` not found in SA vs MOCK up-regulated file.")
}

sa_up_raw <- sa_up_raw |>
  dplyr::mutate(
    Protein_IDs = as.character(.data$Protein_IDs),
    First_ID = sub(";.*", "", .data$Protein_IDs)
  )

if ("Annotated_Gene" %in% colnames(sa_up_raw)) {
  sa_up_raw <- dplyr::select(sa_up_raw, -dplyr::all_of("Annotated_Gene"))
}

sa_up_annot <- sa_up_raw |>
  dplyr::left_join(annot_lookup, by = "First_ID") |>
  dplyr::mutate(
    Annotated_Gene = ifelse(
      is.na(.data$Annotated_Gene) | .data$Annotated_Gene == "",
      .data$First_ID,
      .data$Annotated_Gene
    )
  )

openxlsx::write.xlsx(
  x = sa_up_annot,
  file = file.path(output_dir, "Upregulated_in_SA_SA_vs_MOCK_annotated.xlsx"),
  rowNames = FALSE
)

# ---- [3) Load Top20_Contributing_Comp1.tsv and use it for corrected label overlap] ----
top20_contrib_file <- file.path("Top_Contributing_Proteins", "Top20_Contributing_Comp1.tsv")

if (!file.exists(top20_contrib_file)) {
  stop(
    paste0(
      "File not found: ",
      normalizePath(top20_contrib_file, winslash = "/", mustWork = FALSE)
    )
  )
}

top20_contrib <- readr::read_tsv(
  file = top20_contrib_file,
  show_col_types = FALSE
)

if (!"Gene_names_for_plot" %in% colnames(top20_contrib)) {
  stop("Column `Gene_names_for_plot` not found in Top20_Contributing_Comp1.tsv")
}

# expand Top20 labels from exported table
top20_labels_df <- data.frame(
  Source = "Top20_Comp1",
  Original_Label = top20_contrib$Gene_names_for_plot,
  stringsAsFactors = FALSE
)

top20_labels_df <- do.call(
  rbind,
  lapply(seq_len(nrow(top20_labels_df)), function(i) {
    labs <- split_ids(top20_labels_df$Original_Label[i])
    if (length(labs) == 0) labs <- NA_character_
    data.frame(
      Source = top20_labels_df$Source[i],
      Original_Label = top20_labels_df$Original_Label[i],
      Normalized_Label = labs,
      stringsAsFactors = FALSE
    )
  })
)

top20_labels_df$Normalized_Label <- normalize_labels(top20_labels_df$Normalized_Label)
top20_labels_df <- top20_labels_df |>
  dplyr::filter(!is.na(.data$Normalized_Label), .data$Normalized_Label != "") |>
  dplyr::distinct(.data$Normalized_Label, .keep_all = TRUE)

# build best label for SA rows, then expand
sa_up_annot <- sa_up_annot |>
  dplyr::mutate(
    SA_Best_Label = dplyr::case_when(
      !is.na(.data$Annotated_Gene) & trimws(.data$Annotated_Gene) != "" ~ .data$Annotated_Gene,
      TRUE ~ .data$Protein_IDs
    )
  )

sa_labels_df <- data.frame(
  First_ID = sa_up_annot$First_ID,
  SA_Best_Label = sa_up_annot$SA_Best_Label,
  Protein_IDs = sa_up_annot$Protein_IDs,
  stringsAsFactors = FALSE
)

sa_labels_df <- do.call(
  rbind,
  lapply(seq_len(nrow(sa_labels_df)), function(i) {
    labs <- split_ids(sa_labels_df$SA_Best_Label[i])
    if (length(labs) == 0) labs <- NA_character_
    data.frame(
      First_ID = sa_labels_df$First_ID[i],
      SA_Best_Label = sa_labels_df$SA_Best_Label[i],
      Protein_IDs = sa_labels_df$Protein_IDs[i],
      Normalized_Label = labs,
      stringsAsFactors = FALSE
    )
  })
)

sa_labels_df$Normalized_Label <- normalize_labels(sa_labels_df$Normalized_Label)
sa_labels_df <- sa_labels_df |>
  dplyr::filter(!is.na(.data$Normalized_Label), .data$Normalized_Label != "") |>
  dplyr::distinct(.data$Normalized_Label, .keep_all = TRUE)

# ---- [Corrected overlap: compare normalized labels, not mixed IDs] ----
set1_labels <- sort(unique(top20_labels_df$Normalized_Label))
set2_labels <- sort(unique(sa_labels_df$Normalized_Label))
shared_labels <- sort(intersect(set1_labels, set2_labels))
n_shared <- length(shared_labels)

# ---- [Build corrected shared-protein table WITH ALL columns from data where possible] ----
shared_proteins_df <- data.frame(
  Shared_Label = shared_labels,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(
    top20_labels_df |>
      dplyr::select(
        .data$Normalized_Label,
        Top20_Original_Label = .data$Original_Label
      ),
    by = c("Shared_Label" = "Normalized_Label")
  ) |>
  dplyr::left_join(
    sa_labels_df |>
      dplyr::select(
        .data$Normalized_Label,
        .data$First_ID,
        SA_file_Label = .data$SA_Best_Label,
        SA_file_Protein_IDs = .data$Protein_IDs
      ),
    by = c("Shared_Label" = "Normalized_Label")
  ) |>
  dplyr::left_join(
    top20_annot |>
      dplyr::select(
        .data$First_ID,
        Comp1_Loading = .data$loading,
        Comp1_Annotated_Gene = .data$Annotated_Gene,
        Comp1_Human_GO = .data$Human_GO,
        Comp1_Human_KEGG = .data$Human_KEGG,
        Comp1_Human_Reactome = .data$Human_Reactome
      ) |>
      dplyr::distinct(.data$First_ID, .keep_all = TRUE),
    by = "First_ID"
  ) |>
  dplyr::left_join(
    data_expanded |>
      dplyr::group_by(.data$First_ID) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(all_data_columns),
          ~ collapse_unique(.x)
        ),
        .groups = "drop"
      ),
    by = "First_ID"
  ) |>
  dplyr::mutate(
    Shared_Label_Display = dplyr::case_when(
      !is.na(.data$Shared_Label) & trimws(.data$Shared_Label) != "" ~ .data$Shared_Label,
      !is.na(.data$Comp1_Annotated_Gene) & trimws(.data$Comp1_Annotated_Gene) != "" ~ .data$Comp1_Annotated_Gene,
      TRUE ~ .data$First_ID
    )
  ) |>
  dplyr::arrange(dplyr::desc(.data$Comp1_Loading), .data$Shared_Label_Display)

# ---- [Put useful columns first, then ALL original data columns] ----
front_cols <- c(
  "Shared_Label",
  "Shared_Label_Display",
  "First_ID",
  "Top20_Original_Label",
  "SA_file_Label",
  "SA_file_Protein_IDs",
  "Comp1_Loading",
  "Comp1_Annotated_Gene",
  "Comp1_Human_GO",
  "Comp1_Human_KEGG",
  "Comp1_Human_Reactome"
)

front_cols <- front_cols[front_cols %in% colnames(shared_proteins_df)]
ordered_cols <- unique(c(front_cols, all_data_columns, colnames(shared_proteins_df)))

shared_proteins_df <- shared_proteins_df |>
  dplyr::select(dplyr::all_of(ordered_cols))

# ---- [Write shared proteins outputs] ----
utils::write.table(
  x = shared_proteins_df,
  file = file.path(output_dir, "Shared_Proteins_Top20_Comp1_vs_SA_MOCK_ALL_DATA_COLUMNS.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

openxlsx::write.xlsx(
  x = shared_proteins_df,
  file = file.path(output_dir, "Shared_Proteins_Top20_Comp1_vs_SA_MOCK_ALL_DATA_COLUMNS.xlsx"),
  rowNames = FALSE
)

# ---- [Concise summary export] ----
shared_summary_df <- shared_proteins_df |>
  dplyr::select(
    dplyr::any_of(c(
      "Shared_Label",
      "Shared_Label_Display",
      "First_ID",
      "Protein_IDs",
      "Gene_names",
      "Comp1_Loading",
      "Comp1_Annotated_Gene",
      "Comp1_Human_GO",
      "Comp1_Human_KEGG",
      "Comp1_Human_Reactome"
    ))
  )

utils::write.table(
  x = shared_summary_df,
  file = file.path(output_dir, "Shared_Proteins_Summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# ---- [4) Corrected Venn Diagram] ----
gene_text <- if (n_shared > 0) {
  paste0(seq_along(shared_labels), ". ", shared_labels, collapse = "\n")
} else {
  "No shared proteins"
}

writeLines(
  text = c(
    paste0("Shared proteins (n = ", n_shared, ")"),
    "",
    shared_labels
  ),
  con = file.path(output_dir, "Shared_Proteins_List.txt")
)

vd <- VennDiagram::venn.diagram(
  x = list(
    Top20_Comp1 = set1_labels,
    Signif_SA_vs_MOCK = set2_labels
  ),
  filename = NULL,
  fill = c("#0065BD", "#E41A1C"),
  alpha = 0.35,
  cat.col = c("#0065BD", "#E41A1C"),
  cex = 4.5,
  cat.cex = 3.2,
  cat.pos = c(0, 180),
  cat.dist = c(-0.03, -0.03),
  margin = 0.08
)

# ---- [Preview] ----
grid::grid.newpage()
grid::grid.draw(vd)

grid::grid.text(
  label = paste0("Shared proteins (n = ", n_shared, ")"),
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.16, "npc"),
  just = "center",
  gp = grid::gpar(fontsize = 22, fontface = "bold")
)

grid::grid.text(
  label = gene_text,
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.07, "npc"),
  just = "center",
  gp = grid::gpar(fontsize = 14, fontface = "italic")
)

# ---- [Save PNG] ----
venn_out <- file.path(output_dir, "Venn_Top20_Comp1_vs_SA_shared_improved.png")

grDevices::png(
  filename = venn_out,
  width = 3600,
  height = 4200,
  res = 300
)

grid::grid.newpage()
grid::grid.draw(vd)

grid::grid.text(
  label = paste0("Shared proteins (n = ", n_shared, ")"),
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.16, "npc"),
  just = "center",
  gp = grid::gpar(fontsize = 22, fontface = "bold")
)

grid::grid.text(
  label = gene_text,
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.07, "npc"),
  just = "center",
  gp = grid::gpar(fontsize = 14, fontface = "italic")
)

grDevices::dev.off()

# ---- [Console Report] ----
cat("Annotation and corrected shared-protein analysis complete\n")
cat("─────────────────────────────────────────────────────────\n")
cat("Top20 Comp1 table: ", file.path(output_dir, "Top20_Comp1.txt"), "\n", sep = "")
cat(
  "Top20 Comp1 annotated XLSX: ",
  file.path(output_dir, "Top20_Comp1_annotated.xlsx"),
  "\n",
  sep = ""
)
cat(
  "SA vs MOCK annotated XLSX: ",
  file.path(output_dir, "Upregulated_in_SA_SA_vs_MOCK_annotated.xlsx"),
  "\n",
  sep = ""
)
cat(
  "Shared proteins TSV (all data columns): ",
  file.path(output_dir, "Shared_Proteins_Top20_Comp1_vs_SA_MOCK_ALL_DATA_COLUMNS.tsv"),
  "\n",
  sep = ""
)
cat(
  "Shared proteins XLSX (all data columns): ",
  file.path(output_dir, "Shared_Proteins_Top20_Comp1_vs_SA_MOCK_ALL_DATA_COLUMNS.xlsx"),
  "\n",
  sep = ""
)
cat(
  "Shared proteins summary TSV: ",
  file.path(output_dir, "Shared_Proteins_Summary.tsv"),
  "\n",
  sep = ""
)
cat(
  "Shared protein list TXT: ",
  file.path(output_dir, "Shared_Proteins_List.txt"),
  "\n",
  sep = ""
)
cat("Shared proteins: ", n_shared, "\n", sep = "")
cat("Venn diagram PNG: ", venn_out, "\n", sep = "")
cat("Done.\n")

################################################################################
# STEP 4
# Annotation + enrichment for:
# 1) Shared proteins from volcano workflow
# 2) Top 20 proteins from each sPLS-DA component
# Fully namespaced with package::function usage
#
# FIXES INCLUDED:
# - Species annotation now uses FASTA headers first
# - No longer assumes UniProt-like accessions are Human
# - Better SA detection from FASTA header / protein IDs
# - Adds species count export for transparency
################################################################################

# ---- [Optional Package Checks] ----
required_pkgs_step4 <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "dplyr",
  "openxlsx",
  "AnnotationDbi"
)

missing_pkgs_step4 <- required_pkgs_step4[
  !vapply(required_pkgs_step4, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs_step4) > 0) {
  stop(
    paste0(
      "Please install the following packages before running STEP 4: ",
      paste(missing_pkgs_step4, collapse = ", ")
    )
  )
}

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis")

# ---- [Directories] ----
splsda_dir <- "sPLSDA_2D_Plots"
top20_dir <- "Top_Contributing_Proteins"
step4_dir <- file.path(splsda_dir, "Step4_Annotation_Enrichment")
dir.create(step4_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Input Files] ----
shared_file <- file.path(splsda_dir, "Shared_Proteins_Top20_Comp1_vs_SA_MOCK_ALL_DATA_COLUMNS.xlsx")
top20_combined_file <- file.path(top20_dir, "Top20_Contributing_Proteins_All_Components.tsv")

# ---- [Helper: annotate species] ----
annotate_species <- function(df) {
  if (!"Protein_IDs" %in% colnames(df)) {
    stop("Column `Protein_IDs` not found in input file.")
  }
  
  if (!"Fasta_headers" %in% colnames(df)) {
    df$Fasta_headers <- NA_character_
  }
  
  if (!"Majority_protein_IDs" %in% colnames(df)) {
    df$Majority_protein_IDs <- NA_character_
  }
  
  dplyr::mutate(
    df,
    Protein_IDs = as.character(.data$Protein_IDs),
    Fasta_headers = as.character(.data$Fasta_headers),
    Majority_protein_IDs = as.character(.data$Majority_protein_IDs),
    
    First_ID = sub(";.*", "", .data$Protein_IDs),
    
    Species = dplyr::case_when(
      # ---- explicit prefixed IDs ----
      grepl("^hsa:", .data$First_ID, ignore.case = TRUE) ~ "Human",
      grepl("^sav:|^sau:", .data$First_ID, ignore.case = TRUE) ~ "SA",
      
      # ---- FASTA header-based HUMAN detection ----
      grepl("OX=9606", .data$Fasta_headers, ignore.case = TRUE) ~ "Human",
      grepl("Homo[_ ]sapiens", .data$Fasta_headers, ignore.case = TRUE) ~ "Human",
      grepl("_HUMAN_", .data$Fasta_headers, ignore.case = TRUE) ~ "Human",
      
      # ---- FASTA header-based S. aureus detection ----
      grepl("Staphylococcus[_ ]aureus", .data$Fasta_headers, ignore.case = TRUE) ~ "SA",
      grepl("OX=1280", .data$Fasta_headers, ignore.case = TRUE) ~ "SA",
      grepl("OX=158878", .data$Fasta_headers, ignore.case = TRUE) ~ "SA",
      grepl("_STAAU_|_STAAM_|_STA", .data$Fasta_headers, ignore.case = TRUE) ~ "SA",
      
      # ---- explicit SA-like IDs / names ----
      grepl("^SA", .data$First_ID, ignore.case = TRUE) ~ "SA",
      grepl("^sav:|^sau:", .data$Protein_IDs, ignore.case = TRUE) ~ "SA",
      grepl("^sav:|^sau:", .data$Majority_protein_IDs, ignore.case = TRUE) ~ "SA",
      
      # ---- conservative fallback ----
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
  
  required_cols <- c("Protein_IDs", "First_ID", "Species")
  missing_cols <- required_cols[!required_cols %in% colnames(sig_data)]
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Input data must contain these columns: ",
        paste(required_cols, collapse = ", "),
        ". Missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # ---- [Split IDs by species] ----
  human_ids_raw <- unique(stats::na.omit(sig_data$First_ID[sig_data$Species == "Human"]))
  sa_ids_raw <- unique(stats::na.omit(sig_data$First_ID[sig_data$Species == "SA"]))
  
  # Remove prefixes if present
  human_ids <- gsub("^hsa:", "", human_ids_raw, ignore.case = TRUE)
  sa_ids <- gsub("^sav:|^sau:", "", sa_ids_raw, ignore.case = TRUE)
  
  # ---- [Write input summary] ----
  summary_df <- data.frame(
    Category = c("Total_rows", "Human_IDs_raw", "SA_IDs_raw", "Unknown_IDs_raw"),
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
  
  # ---- [Write species counts] ----
  species_count_df <- as.data.frame(table(sig_data$Species), stringsAsFactors = FALSE)
  colnames(species_count_df) <- c("Species", "Count")
  
  openxlsx::write.xlsx(
    species_count_df,
    file = file.path(out_dir, paste0("Species_counts_", label_suffix, ".xlsx")),
    rowNames = FALSE
  )
  
  # --------------------------------------------------------------------------
  # HUMAN ENRICHMENT
  # --------------------------------------------------------------------------
  valid_uniprot <- character(0)
  human_entrez <- NULL
  
  if (length(human_ids) > 0) {
    valid_uniprot <- human_ids[
      human_ids %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "UNIPROT")
    ]
    
    if (length(valid_uniprot) > 0) {
      human_entrez <- tryCatch(
        {
          clusterProfiler::bitr(
            valid_uniprot,
            fromType = "UNIPROT",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db::org.Hs.eg.db
          )
        },
        error = function(e) {
          cat("Human ID conversion failed for", label_suffix, ":", e$message, "\n")
          NULL
        }
      )
      
      if (!is.null(human_entrez) && nrow(human_entrez) > 0) {
        openxlsx::write.xlsx(
          human_entrez,
          file = file.path(out_dir, paste0("Human_ID_conversion_", label_suffix, ".xlsx")),
          rowNames = FALSE
        )
        
        human_genes <- unique(stats::na.omit(human_entrez$ENTREZID))
        
        # ---- [GO enrichment] ----
        if (length(human_genes) > 0) {
          go_enrich <- tryCatch(
            {
              clusterProfiler::enrichGO(
                gene = human_genes,
                OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE
              )
            },
            error = function(e) {
              cat("GO enrichment failed for", label_suffix, ":", e$message, "\n")
              NULL
            }
          )
          
          if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
            openxlsx::write.xlsx(
              as.data.frame(go_enrich),
              file = file.path(out_dir, paste0("GO_Human_enrichment_", label_suffix, ".xlsx")),
              rowNames = FALSE
            )
          }
          
          # ---- [KEGG enrichment: human] ----
          kegg_human <- tryCatch(
            {
              clusterProfiler::enrichKEGG(
                gene = human_genes,
                organism = "hsa",
                keyType = "kegg",
                pvalueCutoff = 1,
                qvalueCutoff = 1
              )
            },
            error = function(e) {
              cat("Human KEGG enrichment failed for", label_suffix, ":", e$message, "\n")
              NULL
            }
          )
          
          if (!is.null(kegg_human) && nrow(as.data.frame(kegg_human)) > 0) {
            kegg_human_df <- tryCatch(
              {
                as.data.frame(
                  clusterProfiler::setReadable(
                    kegg_human,
                    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                    keyType = "ENTREZID"
                  )
                )
              },
              error = function(e) {
                as.data.frame(kegg_human)
              }
            )
            
            openxlsx::write.xlsx(
              kegg_human_df,
              file = file.path(out_dir, paste0("KEGG_Human_enrichment_", label_suffix, ".xlsx")),
              rowNames = FALSE
            )
          }
          
          # ---- [Reactome enrichment: human] ----
          reactome_human <- tryCatch(
            {
              ReactomePA::enrichPathway(
                gene = human_genes,
                organism = "human",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE
              )
            },
            error = function(e) {
              cat("Reactome enrichment failed for", label_suffix, ":", e$message, "\n")
              NULL
            }
          )
          
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
    sa_kegg <- tryCatch(
      {
        clusterProfiler::enrichKEGG(
          gene = sa_ids,
          organism = "sau",
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      },
      error = function(e) {
        cat("S. aureus KEGG enrichment failed for", label_suffix, ":", e$message, "\n")
        NULL
      }
    )
    
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

# ---- [Process shared proteins from volcano workflow] ----
process_shared_proteins <- function(shared_file, base_out_dir) {
  if (!file.exists(shared_file)) {
    cat("Shared file not found:", shared_file, "\n")
    return(invisible(NULL))
  }
  
  cat("\n==============================\n")
  cat("Processing shared proteins...\n")
  cat("==============================\n")
  
  shared_data <- openxlsx::read.xlsx(shared_file)
  
  if (nrow(shared_data) == 0) {
    cat("Shared protein file is empty.\n")
    return(invisible(NULL))
  }
  
  shared_annot <- annotate_species(shared_data)
  
  shared_out_dir <- file.path(base_out_dir, "annotation_shared")
  dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  annot_file <- file.path(shared_out_dir, "Shared_Proteins_annotated_for_enrichment.xlsx")
  save_annotated_table(shared_annot, annot_file)
  
  run_enrichment_from_df(
    sig_data = shared_annot,
    out_dir = shared_out_dir,
    label_suffix = "shared_Top20_vs_SA_vs_MOCK"
  )
  
  cat("Shared proteins processed.\n")
  invisible(NULL)
}

# ---- [Process top 20 proteins from each component] ----
process_component_top20 <- function(top20_combined_file, top20_dir, base_out_dir) {
  cat("\n===========================================\n")
  cat("Processing top 20 proteins for each component...\n")
  cat("===========================================\n")
  
  top20_data <- NULL
  
  # Priority 1: combined TSV from STEP 2
  if (file.exists(top20_combined_file)) {
    top20_data <- utils::read.delim(
      file = top20_combined_file,
      header = TRUE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    # Fallback: read separate per-component files if present
    comp_files <- list.files(
      top20_dir,
      pattern = "^Top20_Contributing_Comp[0-9]+\\.tsv$",
      full.names = TRUE
    )
    
    if (length(comp_files) > 0) {
      top20_list <- lapply(comp_files, function(f) {
        utils::read.delim(
          file = f,
          header = TRUE,
          sep = "\t",
          quote = "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      })
      top20_data <- dplyr::bind_rows(top20_list)
    }
  }
  
  if (is.null(top20_data)) {
    cat("No STEP 2 top20 input files found.\n")
    cat("Expected either:\n")
    cat("  ", top20_combined_file, "\n", sep = "")
    cat("or per-component TSVs in:\n")
    cat("  ", normalizePath(top20_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
    return(invisible(NULL))
  }
  
  if (!"Component" %in% colnames(top20_data)) {
    stop("Column `Component` not found in top20 data.")
  }
  
  if (!"Protein_IDs" %in% colnames(top20_data)) {
    stop("Column `Protein_IDs` not found in top20 data.")
  }
  
  top20_annot <- annotate_species(top20_data)
  
  # Save combined annotated file
  comp_base_dir <- file.path(base_out_dir, "annotation_top20_components")
  dir.create(comp_base_dir, recursive = TRUE, showWarnings = FALSE)
  
  openxlsx::write.xlsx(
    top20_annot,
    file = file.path(comp_base_dir, "Top20_Contributing_Proteins_All_Components_annotated.xlsx"),
    rowNames = FALSE
  )
  
  components <- unique(top20_annot$Component)
  components <- components[!is.na(components)]
  
  for (comp in components) {
    comp_df <- dplyr::filter(top20_annot, .data$Component == comp)
    
    if (nrow(comp_df) == 0) {
      next
    }
    
    comp_dir <- file.path(comp_base_dir, comp)
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    
    comp_annot_file <- file.path(
      comp_dir,
      paste0("Top20_", comp, "_annotated_for_enrichment.xlsx")
    )
    
    openxlsx::write.xlsx(comp_df, comp_annot_file, rowNames = FALSE)
    
    run_enrichment_from_df(
      sig_data = comp_df,
      out_dir = comp_dir,
      label_suffix = paste0("Top20_", comp)
    )
    
    cat("Processed:", comp, "with", nrow(comp_df), "proteins\n")
  }
  
  invisible(NULL)
}

# ---- [Run shared protein enrichment] ----
process_shared_proteins(
  shared_file = shared_file,
  base_out_dir = step4_dir
)

# ---- [Run component top20 enrichment] ----
process_component_top20(
  top20_combined_file = top20_combined_file,
  top20_dir = top20_dir,
  base_out_dir = step4_dir
)

# ---- [Final report] ----
cat("\n====================================\n")
cat("STEP 4 annotation + enrichment done\n")
cat("====================================\n")
cat("Output folder:\n")
cat(normalizePath(step4_dir, mustWork = FALSE), "\n")
cat("\nGenerated analyses:\n")
cat("1) Shared proteins from volcano workflow\n")
cat("2) Top 20 proteins from each sPLS-DA component\n")
cat("\nSupported enrichment:\n")
cat("- Human: GO, KEGG, Reactome\n")
cat("- S. aureus: KEGG\n")
cat("\nAll STEP 4 analyses completed.\n")

################################################################################
# STEP 5 (FIXED)
# Fisher's exact test with robust gene matching + debug output
################################################################################

# ---- [Libraries] ----
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# ---- [Set Paths] ----
base_dir <- "C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis"

shared_file <- file.path(
  base_dir,
  "sPLSDA_2D_Plots/Step4_Annotation_Enrichment/annotation_shared",
  "Shared_Proteins_annotated_for_enrichment.xlsx"
)

shared_dir <- file.path(
  base_dir,
  "sPLSDA_2D_Plots/Step4_Annotation_Enrichment/annotation_shared"
)

component_annot_dir <- file.path(
  base_dir,
  "sPLSDA_2D_Plots/Step4_Annotation_Enrichment/annotation_top20_components"
)

background_file_txt <- file.path(base_dir, "Adjusted_Normalized_P257_07.txt")

# ---- [Load Background] ----
background_data <- read.table(
  background_file_txt,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

background_data$Gene_names <- ifelse(
  is.na(background_data$Gene_names) | trimws(background_data$Gene_names) == "",
  background_data$Protein_IDs,
  background_data$Gene_names
)

# 🔥 IMPORTANT: normalize gene names
normalize_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == ""] <- NA
  x
}

background_genes <- unique(na.omit(normalize_gene(background_data$Gene_names)))

cat("🧬 Background genes:", length(background_genes), "\n")

# ---- [Helpers] ----
prepare_sig_genes <- function(sig_df) {
  if (!"Annotated_Gene" %in% colnames(sig_df)) {
    sig_df$Annotated_Gene <- NA_character_
  }
  
  if ("Gene_names" %in% colnames(sig_df)) {
    sig_df$Annotated_Gene <- ifelse(
      is.na(sig_df$Annotated_Gene) | trimws(sig_df$Annotated_Gene) == "",
      sig_df$Gene_names,
      sig_df$Annotated_Gene
    )
  }
  
  normalize_gene(sig_df$Annotated_Gene) %>% na.omit() %>% unique()
}

parse_gene_list <- function(x) {
  vals <- strsplit(as.character(x), "[;/,]")[[1]]
  normalize_gene(vals)
}

# ---- [Core Fisher Function] ----
run_fisher_test_set <- function(sig_file, enrichment_dir, output_dir, label_suffix) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!file.exists(sig_file)) {
    cat("❌ Missing:", sig_file, "\n")
    return()
  }
  
  sig_df <- read_xlsx(sig_file)
  sig_genes <- prepare_sig_genes(sig_df)
  
  cat("\n==============================\n")
  cat("📊 DATASET:", label_suffix, "\n")
  cat("Sig genes:", length(sig_genes), "\n")
  
  if (length(sig_genes) == 0) return()
  
  sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
  
  for (source in sources) {
    
    enrich_file <- file.path(
      enrichment_dir,
      paste0(source, "_enrichment_", label_suffix, ".xlsx")
    )
    
    if (!file.exists(enrich_file)) {
      cat("⚠️ Missing:", enrich_file, "\n")
      next
    }
    
    enrich_data <- read_xlsx(enrich_file)
    
    if (!"geneID" %in% colnames(enrich_data)) next
    
    results <- list()
    
    for (i in seq_len(nrow(enrich_data))) {
      
      pathway <- if ("Description" %in% colnames(enrich_data)) {
        enrich_data$Description[i]
      } else {
        paste0("Path_", i)
      }
      
      term_genes <- intersect(parse_gene_list(enrich_data$geneID[i]), background_genes)
      sig_in_pathway <- intersect(term_genes, sig_genes)
      
      a <- length(sig_in_pathway)
      
      # 🔥 skip useless pathways
      if (a == 0) next
      
      b <- length(sig_genes) - a
      c <- length(term_genes) - a
      d <- length(background_genes) - a - b - c
      
      if (any(c(a,b,c,d) < 0)) next
      
      ft <- tryCatch(fisher.test(matrix(c(a,b,c,d),2)), error=function(e) NULL)
      if (is.null(ft)) next
      
      results[[length(results)+1]] <- data.frame(
        Pathway = pathway,
        Sig_In_Pathway = a,
        Genes = paste(sig_in_pathway, collapse=";"),
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
    
    if (length(results) == 0) {
      cat("⚠️ No enrichment for:", source, "\n")
      next
    }
    
    df <- bind_rows(results) %>%
      mutate(p_adj = p.adjust(p_value, method="BH")) %>%
      arrange(p_adj)
    
    cat("✅", source, "rows:", nrow(df), "\n")
    
    write.xlsx(
      df,
      file.path(output_dir, paste0("Fisher_Pathway_", source, "_", label_suffix, ".xlsx")),
      rowNames = FALSE
    )
  }
}

# ---- [RUN SHARED] ----
run_fisher_test_set(
  sig_file = shared_file,
  enrichment_dir = shared_dir,
  output_dir = shared_dir,
  label_suffix = "shared_Top20_vs_SA_vs_MOCK"
)

# ---- [RUN COMPONENTS] ----
comp_dirs <- list.dirs(component_annot_dir, recursive = FALSE)

for (comp_dir in comp_dirs) {
  
  comp <- basename(comp_dir)
  
  run_fisher_test_set(
    sig_file = file.path(comp_dir, paste0("Top20_", comp, "_annotated_for_enrichment.xlsx")),
    enrichment_dir = comp_dir,
    output_dir = file.path(comp_dir, "Fisher_Test"),
    label_suffix = paste0("Top20_", comp)
  )
}

cat("\n✅ FIXED Fisher STEP DONE\n")

################################################################################

################################################################################
# STEP 6
# Enrichment Heatmap Plotting for Fisher's Exact Test Results
# Uses Fisher outputs from STEP 5
# Saves shared plots in annotation_shared
# Saves component plots inside each component/Fisher_Test folder
################################################################################

# ---- [Libraries] ----
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)   # for unit()

# ---- [Set Paths] ----
base_dir <- "C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis"

# Shared annotation directory
shared_dir <- file.path(
  base_dir,
  "sPLSDA_2D_Plots",
  "Step4_Annotation_Enrichment",
  "annotation_shared"
)

# Shared annotated gene file
shared_file <- file.path(
  shared_dir,
  "Shared_Proteins_annotated_for_enrichment.xlsx"
)

# Shared Fisher outputs from STEP 5 are written directly into annotation_shared
shared_fisher_dir <- shared_dir

# Component parent directory
component_annot_dir <- file.path(
  base_dir,
  "sPLSDA_2D_Plots",
  "Step4_Annotation_Enrichment",
  "annotation_top20_components"
)

# ---- [Settings] ----
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
p_adj_cutoff <- 0.05
min_sig_in_pathway <- 3
pathways_per_page <- 20
use_sepsis_filter <- TRUE

sepsis_keywords <- c(
  "sepsis", "infection", "endotoxin", "bacteria", "pathogen", "defense", "acute",
  "innate immune", "acute inflammation response",
  "vesicle", "exosome", "positive regulation of exocytosis", "extracellular vesicle",
  "blood coagulation", "positive regulation of secretion", "blood microparticle",
  "leukocyte chemotaxis", "lytic vacuole membrane", "secretory granule membrane",
  "complement activation",
  "defense response to gram-positive bacterium",
  "positive regulation interleukin-1 production", "macrophage activation",
  "b cell", "neutrophil", "t cell", "adaptive immune"
)

# ---- [Helper: shorten long gene labels] ----
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
  many_ids <- str_count(x, ";") >= 2
  accession_like <- grepl("^[A-Z0-9.;_-]+$", x)
  
  too_long | (many_ids & accession_like)
}

# ---- [Helper: extract significant genes from annotated file] ----
get_sig_genes <- function(annot_file) {
  if (!file.exists(annot_file)) {
    cat("⚠️ Missing annotated file:", annot_file, "\n")
    return(character(0))
  }
  
  annot_df <- read_xlsx(annot_file)
  
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
  
  sig_genes <- unique(trimws(as.character(annot_df$Annotated_Gene)))
  sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]
  
  # Remove problematic labels that break plotting
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

# ---- [Helper: load Fisher result files] ----
load_fisher_results <- function(fisher_dir, label_suffix, sources) {
  combined_all <- list()
  
  for (source in sources) {
    fisher_file <- file.path(
      fisher_dir,
      paste0("Fisher_Pathway_", source, "_", label_suffix, ".xlsx")
    )
    
    if (!file.exists(fisher_file)) {
      cat("⚠️ Missing:", fisher_file, "\n")
      next
    }
    
    fisher_df <- read.xlsx(fisher_file) %>% as.data.frame()
    
    if (nrow(fisher_df) == 0) next
    if (!"Genes" %in% colnames(fisher_df)) next
    if (!"Pathway" %in% colnames(fisher_df)) next
    if (!"p_adj" %in% colnames(fisher_df)) next
    
    sig_count_col <- detect_sig_count_col(fisher_df)
    if (is.na(sig_count_col)) next
    
    fisher_clean <- fisher_df %>%
      filter(
        !is.na(p_adj),
        p_adj < p_adj_cutoff,
        .data[[sig_count_col]] >= min_sig_in_pathway,
        !is.na(Genes),
        Genes != ""
      ) %>%
      mutate(
        Pathway = str_trim(Pathway),
        Source = gsub("_Human", "", source),
        Pathway = paste0(str_wrap(Pathway, 70), " [", Source, "]"),
        log10_p = -log10(p_adj)
      ) %>%
      separate_rows(Genes, sep = "[;,/\\s]+") %>%
      mutate(Gene = str_trim(Genes)) %>%
      filter(!is.na(Gene), Gene != "") %>%
      distinct(Pathway, Gene, log10_p, Source)
    
    combined_all[[source]] <- fisher_clean
  }
  
  bind_rows(combined_all)
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
    cat("⚠️ No Fisher enrichment rows for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  if (use_sepsis_filter) {
    pathway_keep <- unique(
      merged_df$Pathway[
        str_detect(
          tolower(merged_df$Pathway),
          paste(tolower(sepsis_keywords), collapse = "|")
        )
      ]
    )
    
    merged_df <- merged_df %>% filter(Pathway %in% pathway_keep)
    
    if (nrow(merged_df) == 0) {
      cat("⚠️ No sepsis-related pathways found for:", dataset_name, "\n")
      return(invisible(NULL))
    }
  }
  
  # Keep only genes present in Fisher results
  fisher_genes <- sort(unique(merged_df$Gene))
  sig_genes <- intersect(sig_genes, fisher_genes)
  
  if (length(sig_genes) == 0) {
    cat("⚠️ No overlapping genes between annotation and Fisher results for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  # Shortened labels for plotting only
  gene_label_map <- data.frame(
    Gene = sig_genes,
    Gene_Label = clean_gene_label(sig_genes, max_chars = 25),
    stringsAsFactors = FALSE
  )
  
  # Avoid duplicate shortened labels
  gene_label_map$Gene_Label <- make.unique(gene_label_map$Gene_Label)
  
  pathway_levels <- unique(merged_df$Pathway)
  
  full_grid <- expand.grid(
    Gene = sig_genes,
    Pathway = pathway_levels,
    stringsAsFactors = FALSE
  ) %>%
    left_join(merged_df, by = c("Gene", "Pathway")) %>%
    left_join(gene_label_map, by = "Gene")
  
  if (nrow(full_grid) == 0) {
    cat("⚠️ Nothing to plot for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  pathway_ranking <- full_grid %>%
    group_by(Pathway) %>%
    summarize(avg_log10_p = mean(log10_p, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_log10_p))
  
  if (nrow(pathway_ranking) == 0) {
    cat("⚠️ No ranked pathways for:", dataset_name, "\n")
    return(invisible(NULL))
  }
  
  gene_levels <- gene_label_map$Gene_Label[match(sort(sig_genes), gene_label_map$Gene)]
  
  full_grid <- full_grid %>%
    mutate(
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
      filter(Pathway %in% current_paths) %>%
      mutate(Pathway = factor(Pathway, levels = rev(current_paths)))
    
    plot_file <- file.path(
      save_dir,
      paste0(
        gsub("[^A-Za-z0-9_]+", "_", dataset_name),
        "_Heatmap_Page",
        i,
        ".png"
      )
    )
    
    p <- ggplot(chunk_data, aes(x = Gene_Label, y = Pathway, fill = log10_p)) +
      geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
      scale_fill_gradient(
        low = "lightblue",
        high = "red",
        na.value = "black",
        name = "-log10(p.adj)"
      ) +
      labs(
        title = NULL,
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 5) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 32),
        axis.text.y = element_text(size = 32),
        legend.title = element_text(size = 32, margin = margin(b = 15)),
        legend.text = element_text(size = 32),
        legend.key.height = unit(2, "cm"),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(
      filename = plot_file,
      plot = p,
      width = max(30, length(unique(chunk_data$Gene_Label)) * 1),
      height = max(15, length(current_paths) * 0.8),
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    
    cat("✅ Saved:", plot_file, "\n")
  }
  
  invisible(NULL)
}

# ==============================================================================
# STEP 6A: Shared proteins heatmap
# ==============================================================================

shared_sig_genes <- get_sig_genes(shared_file)

if (length(shared_sig_genes) > 0) {
  shared_merged <- load_fisher_results(
    fisher_dir = shared_fisher_dir,
    label_suffix = "shared_Top20_vs_SA_vs_MOCK",
    sources = sources
  )
  
  plot_heatmap_pages(
    merged_df = shared_merged,
    sig_genes = shared_sig_genes,
    dataset_name = "shared_sepsis_EVs_proteins",
    save_dir = shared_fisher_dir,
    use_sepsis_filter = use_sepsis_filter,
    sepsis_keywords = sepsis_keywords,
    pathways_per_page = pathways_per_page
  )
} else {
  cat("⚠️ No shared significant genes found.\n")
}

# ==============================================================================
# STEP 6B: Per-component heatmaps
# ==============================================================================

if (dir.exists(component_annot_dir)) {
  comp_dirs <- list.dirs(component_annot_dir, recursive = FALSE, full.names = TRUE)
  comp_dirs <- comp_dirs[grepl("^Comp[0-9]+$", basename(comp_dirs))]
  
  if (length(comp_dirs) == 0) {
    cat("⚠️ No component directories found in:", component_annot_dir, "\n")
  }
  
  for (comp_dir in comp_dirs) {
    comp_name <- basename(comp_dir)
    
    annot_file <- file.path(
      comp_dir,
      paste0("Top20_", comp_name, "_annotated_for_enrichment.xlsx")
    )
    
    # STEP 5 wrote component Fisher files here:
    fisher_dir <- file.path(comp_dir, "Fisher_Test")
    
    comp_sig_genes <- get_sig_genes(annot_file)
    
    if (length(comp_sig_genes) == 0) {
      cat("⚠️ No significant genes found for", comp_name, "\n")
      next
    }
    
    if (!dir.exists(fisher_dir)) {
      cat("⚠️ Fisher_Test folder not found for", comp_name, ":", fisher_dir, "\n")
      next
    }
    
    comp_merged <- load_fisher_results(
      fisher_dir = fisher_dir,
      label_suffix = paste0("Top20_", comp_name),
      sources = sources
    )
    
    plot_heatmap_pages(
      merged_df = comp_merged,
      sig_genes = comp_sig_genes,
      dataset_name = paste0("Top20_", comp_name),
      save_dir = fisher_dir,
      use_sepsis_filter = use_sepsis_filter,
      sepsis_keywords = sepsis_keywords,
      pathways_per_page = pathways_per_page
    )
  }
} else {
  cat("⚠️ Component annotation directory not found:", component_annot_dir, "\n")
}

# ---- [Final Report] ----
cat("\n=================================================\n")
cat("STEP 6 enrichment heatmap plotting completed\n")
cat("=================================================\n")
cat("Plots saved directly inside each Fisher output folder.\n")
cat("Shared plots saved inside annotation_shared.\n")
cat("Component plots saved inside each component/Fisher_Test folder.\n")
cat("Long accession-style labels were filtered/shortened for plotting.\n")
cat("\n✅ STEP 6 done.\n")

################################################################################

# =============================================================================
# Full volcano plots with top 20 sPLS-DA Comp1 proteins highlighted
# Saves outputs inside: sPLSDA_2D_Plots
# Uses paired t-test
# =============================================================================

# ---- [Libraries] ----
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/adj_p_value_P257_07_new_analysis")

# ---- [Output Directory] ----
main_output_dir <- "sPLSDA_2D_Plots"
dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Checks] ----
if (!exists("splsda_model")) {
  stop("Object `splsda_model` not found. Please run STEP 1 first in the same R session.")
}
if (!exists("data")) {
  stop("Object `data` not found. Please run STEP 1 first in the same R session.")
}
if (is.null(splsda_model$loadings$X) || ncol(splsda_model$loadings$X) < 1) {
  stop("`splsda_model$loadings$X` is missing or has no components.")
}

# ---- [Load Annotation File] ----
Annotation_Human_SA <- read.table(
  "Annotations(Human+all S. aureus).txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))
stopifnot(all(c("Entry", "Gene_Names__primary_") %in% colnames(Annotation_Human_SA)))

# ---- [Gene Cleanup] ----
data$Gene_names <- ifelse(
  is.na(data$Gene_names) | trimws(data$Gene_names) == "",
  data$Protein_IDs,
  data$Gene_names
)
data$Gene_names <- make.unique(as.character(data$Gene_names))

# ---- [Annotation Map] ----
annotation_map <- Annotation_Human_SA %>%
  dplyr::select(Entry, Gene_Names__primary_) %>%
  dplyr::rename(
    First_ID = Entry,
    Annotated_Gene = Gene_Names__primary_
  ) %>%
  dplyr::mutate(
    First_ID = trimws(as.character(First_ID)),
    Annotated_Gene = trimws(as.character(Annotated_Gene))
  ) %>%
  dplyr::distinct(First_ID, .keep_all = TRUE)

# ---- [Map sPLS-DA feature IDs back to original data] ----
data_export <- data
data_export$Protein_IDs <- as.character(data_export$Protein_IDs)
data_export$Gene_names <- as.character(data_export$Gene_names)

data_export$Gene_names <- ifelse(
  is.na(data_export$Gene_names) | trimws(data_export$Gene_names) == "",
  data_export$Protein_IDs,
  data_export$Gene_names
)

data_export$Gene_names <- make.unique(as.character(data_export$Gene_names))
data_export$Feature_ID_for_join <- make.names(data_export$Gene_names, unique = TRUE)

# ---- [Extract top 20 proteins from Comp1] ----
top_n_splsda <- 20
comp_to_use <- 1

loadings_vec <- splsda_model$loadings$X[, comp_to_use]
loadings_vec <- loadings_vec[!is.na(loadings_vec)]

ord <- order(abs(loadings_vec), decreasing = TRUE)
top_features <- loadings_vec[ord][seq_len(min(top_n_splsda, length(loadings_vec)))]

top20_df <- data.frame(
  Feature_ID_for_join = names(top_features),
  Loading = as.numeric(top_features),
  Abs_Loading = abs(as.numeric(top_features)),
  Rank = seq_along(top_features),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(data_export, by = "Feature_ID_for_join") %>%
  dplyr::mutate(
    First_ID = sub(";.*", "", Protein_IDs)
  ) %>%
  dplyr::left_join(annotation_map, by = "First_ID") %>%
  dplyr::mutate(
    Annotated_Gene = ifelse(
      is.na(Annotated_Gene) | trimws(Annotated_Gene) == "",
      First_ID,
      Annotated_Gene
    ),
    Label_for_plot = dplyr::case_when(
      !is.na(Gene_names) & trimws(Gene_names) != "" & Gene_names != Protein_IDs ~ Gene_names,
      !is.na(Annotated_Gene) & trimws(Annotated_Gene) != "" ~ Annotated_Gene,
      TRUE ~ Protein_IDs
    )
  )

write.table(
  top20_df,
  file = file.path(main_output_dir, "Top20_sPLSDA_Comp1_highlighted_in_volcano.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

top20_feature_ids <- unique(top20_df$Feature_ID_for_join)

# ---- [Thresholds] ----
pval_threshold <- 0.05
fc_threshold <- 0.3010  # log10(2)

# ---- [Manual Comparisons] ----
comparisons <- list(
  list(
    name = "SA_vs_MOCK",
    group1_name = "SA",
    group2_name = "MOCK",
    group1_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6"),
    group2_cols = c("MOCK_1", "MOCK_2", "MOCK_3", "MOCK_4", "MOCK_5", "MOCK_6")
  ),
  list(
    name = "PipTazo_vs_SA",
    group1_name = "PipTazo",
    group2_name = "SA",
    group1_cols = c("PipTazo_1", "PipTazo_2", "PipTazo_3", "PipTazo_4", "PipTazo_5", "PipTazo_6"),
    group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6")
  ),
  list(
    name = "Vancomycin_vs_SA",
    group1_name = "Vancomycin",
    group2_name = "SA",
    group1_cols = c("Vancomycin_1", "Vancomycin_2", "Vancomycin_3", "Vancomycin_4", "Vancomycin_5", "Vancomycin_6"),
    group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6")
  ),
  list(
    name = "Moxifloxacin_vs_SA",
    group1_name = "Moxifloxacin",
    group2_name = "SA",
    group1_cols = c("Moxifloxacin_1", "Moxifloxacin_2", "Moxifloxacin_3", "Moxifloxacin_4", "Moxifloxacin_5", "Moxifloxacin_6"),
    group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6")
  )
)

# ---- [Process Each Comparison] ----
for (comp in comparisons) {
  
  comp_name <- comp$name
  group1 <- comp$group1_name
  group2 <- comp$group2_name
  g1_cols <- comp$group1_cols
  g2_cols <- comp$group2_cols
  
  label_up <- paste("Higher in", group1)
  label_down <- paste("Higher in", group2)
  
  df <- data.frame(
    Gene = data$Gene_names,
    Protein_IDs = data$Protein_IDs,
    Feature_ID_for_join = data_export$Feature_ID_for_join,
    data[, c(g1_cols, g2_cols)],
    stringsAsFactors = FALSE
  )
  
  df$mean_diff <- NA_real_
  df$p_value <- NA_real_
  
  for (i in seq_len(nrow(df))) {
    vals1 <- as.numeric(unlist(df[i, g1_cols]))
    vals2 <- as.numeric(unlist(df[i, g2_cols]))
    
    if (all(is.na(vals1)) || all(is.na(vals2))) next
    if (length(vals1) != length(vals2)) next
    
    test_result <- tryCatch(
      t.test(vals1, vals2, paired = TRUE),
      error = function(e) NULL
    )
    
    if (!is.null(test_result)) {
      df$mean_diff[i] <- mean(vals1, na.rm = TRUE) - mean(vals2, na.rm = TRUE)
      df$p_value[i] <- test_result$p.value
    }
  }
  
  # ---- [BH correction] ----
  df$p_value_adj <- p.adjust(df$p_value, method = "BH")
  df$negLog10P <- -log10(df$p_value_adj)
  
  # ---- [Annotation] ----
  df <- df %>%
    dplyr::mutate(
      First_ID = sub(";.*", "", Protein_IDs)
    ) %>%
    dplyr::left_join(annotation_map, by = "First_ID")
  
  df$Annotated_Gene[is.na(df$Annotated_Gene) | trimws(df$Annotated_Gene) == ""] <-
    df$First_ID[is.na(df$Annotated_Gene) | trimws(df$Annotated_Gene) == ""]
  
  df <- df %>%
    dplyr::mutate(
      Label_for_plot = dplyr::case_when(
        !is.na(Gene) & trimws(Gene) != "" & Gene != Protein_IDs ~ Gene,
        !is.na(Annotated_Gene) & trimws(Annotated_Gene) != "" ~ Annotated_Gene,
        TRUE ~ First_ID
      )
    )
  
  # ---- [Volcano grouping] ----
  df$group <- dplyr::case_when(
    df$mean_diff > fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_up,
    df$mean_diff < -fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_down,
    TRUE ~ "Non-significant"
  )
  
  # ---- [Highlight top20 proteins] ----
  df$highlight_top20 <- df$Feature_ID_for_join %in% top20_feature_ids
  
  df$plot_group <- dplyr::case_when(
    df$highlight_top20 & df$group == label_up ~ label_up,
    df$highlight_top20 & df$group == label_down ~ label_down,
    df$highlight_top20 & df$group == "Non-significant" ~ "Non-significant",
    TRUE ~ "Background"
  )
  
  df$plot_group <- factor(
    df$plot_group,
    levels = c(label_up, label_down, "Non-significant", "Background")
  )
  
  # ---- [Labels only for highlighted top20 significant proteins] ----
  label_df <- df %>%
    dplyr::filter(
      highlight_top20,
      group != "Non-significant",
      !is.na(Label_for_plot)
    ) %>%
    dplyr::arrange(dplyr::desc(negLog10P), dplyr::desc(abs(mean_diff)))
  
  # ---- [Save tables] ----
  write.table(
    df,
    file = file.path(main_output_dir, paste0("All_proteins_", comp_name, "_with_Top20_highlight.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    df[df$highlight_top20, ],
    file = file.path(main_output_dir, paste0("Top20_highlighted_proteins_", comp_name, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # ---- [Color map: same style] ----
  color_map <- stats::setNames(
    object = c("red", "blue", "black", "grey70"),
    nm = c(label_up, label_down, "Non-significant", "Background")
  )
  
  # ---- [Volcano Plot] ----
  volcano_plot <- ggplot(df, aes(x = mean_diff, y = negLog10P)) +
    geom_point(aes(color = plot_group), alpha = 0.65, size = 4.2) +
    scale_color_manual(values = color_map, drop = FALSE) +
    geom_hline(
      yintercept = -log10(pval_threshold),
      linetype = "dashed",
      color = "black"
    ) +
    geom_vline(
      xintercept = c(-fc_threshold, fc_threshold),
      linetype = "dashed",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Label_for_plot, color = group),
      size = 7,
      max.overlaps = 100,
      box.padding = 0.5,
      point.padding = 0.3,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title = paste0("Volcano Plot: ", comp_name),
      subtitle = "All proteins, top 20 sPLS-DA Comp1 proteins highlighted",
      x = "Mean difference",
      y = "-Log10 adjusted p-value"
    ) +
    theme_minimal(base_size = 25) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    guides(color = guide_legend(override.aes = list(size = 6)))
  
  print(volcano_plot)
  
  ggsave(
    filename = file.path(main_output_dir, paste0("Volcano_", comp_name, "_Top20_highlight.png")),
    plot = volcano_plot,
    dpi = 300,
    width = 12,
    height = 12,
    bg = "white"
  )
  
  cat("✅ Saved volcano for ", comp_name, "\n", sep = "")
}

cat("🎉 All volcano plots completed.\n")
cat("Saved in: ", normalizePath(main_output_dir), "\n", sep = "")
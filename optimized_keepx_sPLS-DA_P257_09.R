# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots_P257_09"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Load Expression Data] ----
data <- read.table("08.12.25_NEW_Adjusted_Normalized_sepsis_P257_09.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# ---- [Gene Names Fallback] ----
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "", 
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))
rownames(data) <- data$Gene_names

# ---- [Define Group Columns] ----
group_cols <- list(
  HC = grep("^HC_[1-6]$", colnames(data), value = TRUE),
  BCN = grep("^BCN_[1-6]$", colnames(data), value = TRUE),
  BCP = grep("^BCP_[1-6]$", colnames(data), value = TRUE)
)

# ---- [Subset Groups] ----
selected_groups <- names(group_cols)
selected_cols <- unlist(group_cols[selected_groups])
X_data <- data[, selected_cols]
X <- as.data.frame(t(X_data))

# ---- [Assign Group Labels] ----
sample_names <- rownames(X)
Y <- case_when(
  grepl("^HC", sample_names) ~ "HC",
  grepl("^BCN", sample_names) ~ "BCN",
  grepl("^BCP", sample_names) ~ "BCP",
  TRUE ~ NA_character_
)
valid_idx <- which(!is.na(Y))
X <- X[valid_idx, ]
Y <- factor(Y[valid_idx], levels = c("HC", "BCN", "BCP"))

# ---- [Diagnostic Check for Groups] ----
cat("📌 Group Summary (Y):\n")
print(table(Y))
if (length(unique(Y)) < 2) {
  stop("❌ ERROR: 'Y' must contain at least two unique groups.")
}

# ---- [Group Colors] ----
group_colors <- c(
  "HC" = "#E41A1C",
  "BCN" = "#377EB8",
  "BCP" = "#4DAF4A"
)

# ---- [Theme Function] ----
custom_theme_splsda <- function(base_size = 24) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.position = "bottom"
    )
}

# ---- [Tuning sPLS-DA] ----
set.seed(1)
optimal_ncomp <- 2
list_keepX <- c(25, 50, 100)

tune_result <- tune.splsda(
  X, Y,
  ncomp = optimal_ncomp,
  validation = "Mfold",
  folds = 5,
  dist = "centroids.dist",
  measure = "BER",
  test.keepX = list_keepX,
  nrepeat = 10,
  progressBar = TRUE
)

optimal_keepX <- tune_result$choice.keepX[1:optimal_ncomp]
cat("\n✅ Optimal keepX values per component:\n")
print(optimal_keepX)

# ---- [keepX Barplot] ----
keepx_df <- data.frame(Component = factor(seq_along(optimal_keepX)), keepX = optimal_keepX)

gg_keepx <- ggplot(keepx_df, aes(x = Component, y = keepX)) +
  geom_bar(stat = "identity", fill = "#377EB8", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 5) +
  labs(title = "Optimal keepX per Component", x = "Component", y = "Variables Selected") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA))

print(gg_keepx)
ggsave(file.path(output_dir, "keepX_per_component_barplot.png"),
       plot = gg_keepx, dpi = 300, width = 6, height = 5, bg = "white")

# ---- [Fit Final Model] ----
splsda_model <- splsda(X, Y, ncomp = optimal_ncomp, keepX = optimal_keepX)

# ---- [Model Performance] ----
perf_result <- perf(splsda_model, validation = "Mfold", folds = 5, nrepeat = 10,
                    dist = "centroids.dist", progressBar = TRUE)
print(perf_result$error.rate)

# ---- [Explained Variance] ----
expl_var <- apply(splsda_model$variates$X^2, 2, sum) / sum(splsda_model$X^2)

# ---- [Manual Ellipse Function] ----
desired_ellipse_level <- 0.95
compute_ellipse <- function(mean, cov, level = desired_ellipse_level, npoints = 100) {
  angles <- seq(0, 2 * pi, length.out = npoints)
  radius <- sqrt(qchisq(level, df = 2))
  eig <- eigen(cov)
  axes <- radius * t(eig$vectors %*% diag(sqrt(eig$values)))
  ellipse <- t(axes %*% rbind(cos(angles), sin(angles))) + matrix(rep(mean, each = npoints), ncol = 2, byrow = FALSE)
  df <- as.data.frame(ellipse)
  colnames(df) <- c("comp1", "comp2")
  return(df)
}

# ---- [Plot All Component Pairs] ----
for (i in 1:(optimal_ncomp - 1)) {
  for (j in (i + 1):optimal_ncomp) {
    
    plot_data <- data.frame(
      comp1 = splsda_model$variates$X[, i],
      comp2 = splsda_model$variates$X[, j],
      Group = Y,
      Sample = rownames(X)
    )
    
    group_centroids <- plot_data %>%
      group_by(Group) %>%
      summarise(comp1 = mean(comp1), comp2 = mean(comp2), count = n(), .groups = "drop")
    
    ellipse_data <- plot_data %>%
      group_by(Group) %>%
      do({
        group_data <- select(., comp1, comp2)
        ell <- compute_ellipse(colMeans(group_data), cov(group_data))
        ell$Group <- unique(.$Group)
        ell
      }) %>% ungroup()
    
    x_lab <- paste0("Component ", i, " (", round(expl_var[i] * 100, 1), "%)")
    y_lab <- paste0("Component ", j, " (", round(expl_var[j] * 100, 1), "%)")
    
    p <- ggplot(plot_data, aes(x = comp1, y = comp2, color = Group, fill = Group)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_polygon(data = ellipse_data, aes(group = Group), alpha = 0.2, color = NA) +
      geom_path(data = ellipse_data, aes(group = Group), linewidth = 1) +
      geom_text_repel(data = group_centroids,
                      aes(label = paste0(Group, "\n(n=", count, ")")),
                      color = "black", size = 10, fontface = "bold",
                      max.overlaps = 100, box.padding = 0.6, point.padding = 0.6) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      labs(
        title = paste0("sPLS-DA"),
        x = x_lab, y = y_lab
      ) +
      custom_theme_splsda()
    
    print(p)
    
    ggsave(file.path(output_dir, paste0("sPLS-DA_P257_09_Comp", i, "_vs_Comp", j, ".png")),
           plot = p, dpi = 300, width = 10, height = 10, bg = "white")
  }
}

# ---- [Summary Report] ----
cat("\n📊 sPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("✔ Number of Components Used:", optimal_ncomp, "\n")
cat("✔ keepX per Component:", paste(optimal_keepX, collapse = ", "), "\n")
cat("✔ Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")
cat("✔ Explained Variance (First 2 components):", 
    paste0(round(expl_var[1:2] * 100, 1), collapse = "%, "), "%\n")
cat("✔ Output folder:", output_dir, "\n")


################################################################################
################################################################################
################################################################################

# ---- [Load Libraries] ----
library(ggplot2)
library(stringr)
library(dplyr)
library(readr)

# ---- [Set Parameters] ----
top_n <- 20
n_components_to_plot <- 2
annotation_file <- "Annotations(Human+all S. aureus).txt"
output_dir <- "Top_Contributing_Proteins_Plots_P257_09"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Custom Color Palette for Each Component] ----
component_colors <- c(
  "Comp1" = "#F8766D",
  "Comp2" = "#7CAE00"
)

# ---- [Read Annotation File] ----
Annotation_Human_SA <- read.table(annotation_file,
                                  header = TRUE,
                                  sep = "\t",
                                  quote = "",
                                  fill = TRUE,
                                  comment.char = "")
colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", ".", colnames(Annotation_Human_SA))

annotation_map <- Annotation_Human_SA %>%
  select(Entry, Gene.Names..primary., Organism) %>%
  rename(
    First_ID = Entry,
    Annotated_Gene = Gene.Names..primary.,
    Organism_Source = Organism
  )

# ---- [Extract Top Features from sPLS-DA Model] ----
loading_list <- lapply(1:n_components_to_plot, function(comp) {
  loadings <- splsda_model$loadings$X[, comp]
  top_features <- sort(abs(loadings), decreasing = TRUE)[1:top_n]
  genes <- names(top_features)
  
  data.frame(
    Protein_IDs = genes,
    Loading = loadings[genes],
    Component = paste0("Comp", comp),
    stringsAsFactors = FALSE
  )
})
top_loadings_df <- do.call(rbind, loading_list)

# ---- [Annotate with S. aureus info] ----
top_loadings_df <- top_loadings_df %>%
  mutate(First_ID = sub(";.*", "", Protein_IDs)) %>%
  left_join(annotation_map, by = "First_ID") %>%
  mutate(
    Display_Label = case_when(
      !is.na(Organism_Source) & grepl("aureus", Organism_Source, ignore.case = TRUE) ~
        paste0("S. aureus: ", Annotated_Gene),
      TRUE ~ str_trunc(Protein_IDs, 30)
    )
  )

# ---- [Factor Label Order] ----
top_loadings_df <- top_loadings_df %>%
  group_by(Component) %>%
  mutate(Display_Label = factor(Display_Label, levels = rev(unique(Display_Label)))) %>%
  ungroup()

# ---- [Dynamic Plot Settings] ----
base_size <- 14
plot_height <- 5 + n_components_to_plot * 3.5

# ---- [Save Individual Plots by Component] ----
unique_components <- unique(top_loadings_df$Component)

for (comp in unique_components) {
  comp_df <- top_loadings_df %>% filter(Component == comp)
  comp_color <- component_colors[comp]
  
  gg <- ggplot(comp_df, aes(x = Display_Label, y = Loading)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.9, fill = comp_color) +
    coord_flip() +
    labs(
      title = paste0("Top ", top_n, " Contributing Proteins - ", comp),
      x = NULL,
      y = "Loading Value"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.y = element_text(size = base_size - 1, margin = margin(r = 10)),
      axis.text.x = element_text(size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 4),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "none"
    )
  
  print(gg)
  
  ggsave(filename = file.path(output_dir, paste0("Top_Contributing_", comp, ".png")),
         plot = gg, dpi = 300, width = 10, height = 6, bg = "white")
}

# ---- [Save Combined Plot] ----
gg_combined <- ggplot(top_loadings_df, aes(x = Display_Label, y = Loading, fill = Component)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
  coord_flip() +
  facet_wrap(~Component, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_fill_manual(values = component_colors) +
  labs(
    title = paste0("Top ", top_n, " Contributing Proteins per Component"),
    x = NULL,
    y = "Loading Value"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    strip.text = element_text(face = "bold", size = base_size + 2),
    axis.text.y = element_text(size = base_size - 1, margin = margin(r = 10)),
    axis.text.x = element_text(size = base_size - 1),
    axis.title = element_text(face = "bold", size = base_size),
    plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing.y = unit(3.5, "lines"),
    legend.position = "none",
    strip.placement = "outside"
  )

print(gg_combined)

ggsave(filename = file.path(output_dir, "Top_Contributing_All_Components.png"),
       plot = gg_combined, dpi = 300, width = 13, height = plot_height, bg = "white")

cat(paste0("\n✅ Saved top ", top_n, " proteins for 2 components. Check folder: ", output_dir, "\n"))

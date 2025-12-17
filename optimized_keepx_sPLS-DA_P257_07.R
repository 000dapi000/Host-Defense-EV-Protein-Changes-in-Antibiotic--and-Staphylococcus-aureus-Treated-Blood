# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Load Expression Data] ----
data <- read.table("Adjusted_Normalized_P257_07.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "", 
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))
rownames(data) <- data$Gene_names

# ---- [Define Groups] ----
group_cols <- list(
  MOCK = grep("^MOCK_[1-6]$", colnames(data), value = TRUE),
  SA = grep("^SA_[1-6]$", colnames(data), value = TRUE),
  Pip_Tazo = grep("^Pip_Tazo_[1-6]$", colnames(data), value = TRUE),
  Vancomycin = grep("^Vancomycin_[1-6]$", colnames(data), value = TRUE),
  Moxifloxacin = grep("^Moxifloxacin_[1-6]$", colnames(data), value = TRUE)
)

selected_groups <- names(group_cols)
X_data <- data[, unlist(group_cols)]
X <- as.data.frame(t(X_data))

sample_names <- rownames(X)
Y <- case_when(
  grepl("^MOCK", sample_names) ~ "MOCK",
  grepl("^SA", sample_names) ~ "SA",
  grepl("^Pip_Tazo", sample_names) ~ "Pip_Tazo",
  grepl("^Vancomycin", sample_names) ~ "Vancomycin",
  grepl("^Moxifloxacin", sample_names) ~ "Moxifloxacin",
  TRUE ~ NA_character_
)

valid_idx <- which(!is.na(Y))
X <- X[valid_idx, ]
Y <- factor(Y[valid_idx], levels = selected_groups)

group_colors <- c(
  "MOCK" = "#E41A1C", "SA" = "#377EB8", "Pip_Tazo" = "#4DAF4A",
  "Vancomycin" = "#984EA3", "Moxifloxacin" = "#FF7F00"
)

custom_theme_splsda <- function(base_size = 14) {
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
optimal_ncomp <- 4
list_keepX <- c(25, 50, 100)

cat("🔍 Tuning sPLS-DA...\n")
tune_result <- tune.splsda(X, Y, ncomp = optimal_ncomp,
                           validation = "Mfold", folds = 5,
                           dist = "centroids.dist", measure = "BER",
                           test.keepX = list_keepX, nrepeat = 10, progressBar = TRUE)

optimal_keepX <- tune_result$choice.keepX[1:optimal_ncomp]
print(optimal_keepX)

# ---- [keepX Barplot] ----
keepx_df <- data.frame(
  Component = factor(paste0("Comp", seq_along(optimal_keepX))),
  keepX = optimal_keepX
)

gg_keepx <- ggplot(keepx_df, aes(x = Component, y = keepX)) +
  geom_bar(stat = "identity", fill = "#377EB8", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 5) +
  labs(title = "Optimal keepX per Component", x = "Component", y = "Variables Selected") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA))

print(gg_keepx)

ggsave(file.path(output_dir, "keepX_per_component_P257_07.png"),
       plot = gg_keepx, dpi = 300, width = 6, height = 5, bg = "white")

# ---- [Fit Final Model] ----
splsda_model <- splsda(X, Y, ncomp = optimal_ncomp, keepX = optimal_keepX)

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
                      color = "black", size = 5, fontface = "bold",
                      max.overlaps = 100, box.padding = 0.6, point.padding = 0.6) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      labs(
        title = paste0("sPLS-DA: P257_07 (Comp ", i, " vs ", j, ")"),
        x = x_lab, y = y_lab
      ) +
      custom_theme_splsda()
    
    print(p)
    
    ggsave(file.path(output_dir, paste0("sPLS-DA_P257_07_Comp", i, "_vs_Comp", j, ".png")),
           plot = p, dpi = 300, width = 10, height = 10, bg = "white")
  }
}

cat("\n✅ All sPLS-DA 2D plots saved in folder:", output_dir, "\n")

# ---- [Summary Report] ----
cat("\n📊 sPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("✔ Number of Components Used:", optimal_ncomp, "\n")
cat("✔ keepX per Component:", paste(optimal_keepX, collapse = ", "), "\n")
cat("✔ Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")

# Calculate explained variance for first two components or fewer if not available
n_to_report <- min(2, length(expl_var))
cat("✔ Explained Variance (First", n_to_report, "components):",
    paste0(round(expl_var[1:n_to_report] * 100, 1), collapse = "%, "), "%\n")
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
annotation_file <- "Annotations(Human+all S. aureus).txt"
n_components_to_plot <- 4
output_dir <- "Top_Contributing_Proteins_Plots"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Custom Color Palette for Each Component] ----
component_colors <- c(
  "Comp1" = "#F8766D",
  "Comp2" = "#7CAE00",
  "Comp3" = "#00BFC4",
  "Comp4" = "#C77CFF",
  "Comp5" = "#FF61C3",
  "Comp6" = "#00BA38"
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
  
  gg <- ggplot(comp_df, aes(x = Display_Label, y = Loading, fill = Component)) +
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

# ---- [Save Combined Plot if More Than 1 Component] ----
if (n_components_to_plot > 1) {
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
         plot = gg_combined, dpi = 300, width = 10, height = plot_height, bg = "white")
}

cat(paste0("\n✅ Saved top ", top_n, " proteins across ", n_components_to_plot, 
           " components. Check folder: ", output_dir, "\n"))


################################################################################
################################################################################
################################################################################

# ---- [Load Libraries] ----
library(mixOmics)      # for splsda_model
library(dplyr)         # data-wrangling
library(readr)         # fast TSV reading
library(openxlsx)      # write.xlsx
library(VennDiagram)   # venn.diagram
library(grid)          # grid.newpage(), grid.draw(), grid.text()

# ---- [Working Directory & Output Folder] ----
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")
output_dir <- "sPLSDA_2D_Plots"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Load & Sanitize Global Annotation Table] ----
Annotation_Human_SA <- read.table(
  "Annotations(Human+all S. aureus).txt",
  header     = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = ""
)
colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))

# ---- [Build a Clean Annotation Lookup] ----
annot_lookup <- Annotation_Human_SA %>%
  select(
    Entry,
    Gene_Names__primary_,
    Gene_Ontology__GO_,
    KEGG,
    Reactome
  ) %>%
  rename(
    First_ID       = Entry,
    Annotated_Gene = Gene_Names__primary_,
    Human_GO       = Gene_Ontology__GO_,
    Human_KEGG     = KEGG,
    Human_Reactome = Reactome
  )

# ---- [1) Annotate Top-20 Comp1 Contributors] ----
df_ld <- data.frame(
  raw_ID  = names(splsda_model$loadings$X[,1]),
  loading = abs(splsda_model$loadings$X[,1]),
  stringsAsFactors = FALSE
) %>%
  mutate(First_ID = sub(";.*", "", raw_ID))

top20_df <- df_ld %>%
  group_by(First_ID) %>%
  slice_max(order_by = loading, n = 1) %>%
  ungroup() %>%
  arrange(desc(loading)) %>%
  slice_head(n = 20)

top20_annot <- top20_df %>%
  left_join(annot_lookup, by = "First_ID") %>%
  mutate(
    Annotated_Gene = ifelse(
      is.na(Annotated_Gene) | Annotated_Gene == "",
      First_ID,
      Annotated_Gene
    )
  )

write.table(
  top20_df,
  file      = file.path(output_dir, "Top20_Comp1.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
write.xlsx(
  top20_annot,
  file     = file.path(output_dir, "Top20_Comp1_annotated.xlsx"),
  rowNames = FALSE
)

# ---- [2) Annotate SA vs MOCK Up-regulated List] ----
sa_up_file <- file.path("SA_vs_MOCK", "Upregulated_in_SA_SA_vs_MOCK.txt")
sa_up_raw  <- read_tsv(sa_up_file, show_col_types = FALSE) %>%
  mutate(First_ID = sub(";.*", "", Protein_IDs))

# Drop any pre-existing Annotated_Gene column to avoid conflicts
if ("Annotated_Gene" %in% colnames(sa_up_raw)) {
  sa_up_raw <- sa_up_raw %>% select(-Annotated_Gene)
}

sa_up_annot <- sa_up_raw %>%
  left_join(annot_lookup, by = "First_ID") %>%
  mutate(
    Annotated_Gene = ifelse(
      is.na(Annotated_Gene) | Annotated_Gene == "",
      First_ID,
      Annotated_Gene
    )
  )

write.xlsx(
  sa_up_annot,
  file     = file.path(output_dir, "Upregulated_in_SA_SA_vs_MOCK_annotated.xlsx"),
  rowNames = FALSE
)

# ---- [3) Draw & Save Venn Diagram] ----
set1 <- top20_annot$Annotated_Gene
set2 <- unique(sa_up_annot$Annotated_Gene)

shared_proteins <- intersect(set1, set2)
legend_str_lines <- c(
  "Shared proteins between Top20_Comp1 and SA_vs_MOCK:", 
  paste(shared_proteins, collapse = ", ")
)

vd <- venn.diagram(
  x        = list(Top20_Comp1       = set1,
                  Signif_SA_vs_MOCK = set2),
  filename = NULL,
  fill     = c("#377EB8", "#E41A1C"),
  alpha    = 0.5,
  cat.col  = c("#377EB8", "#E41A1C"),
  cex      = 3.5,   # Venn count size
  cat.cex  = 3.5,   # Venn category label size
  cat.pos  = c(-20, 20)
)

# --- Preview in R ---
grid.newpage()
grid.text(legend_str_lines[1], x = unit(0.5, "npc"), y = unit(0.96, "npc"),
          just = "center", gp = gpar(fontsize = 18, fontface = "bold"))
grid.draw(vd)
grid.text(legend_str_lines[2], x = unit(0.5, "npc"), y = unit(0.04, "npc"),
          just = "center", gp = gpar(fontsize = 16, fontface = "italic"))

# --- Save high-res PNG ---
venn_out <- file.path(output_dir, "Venn_Top20_Comp1_vs_SA_shared_horizontal.png")
png(venn_out, width = 3600, height = 3600, res = 300)
grid.newpage()
grid.text(legend_str_lines[1], x = unit(0.5, "npc"), y = unit(0.96, "npc"),
          just = "center", gp = gpar(fontsize = 22, fontface = "bold"))
grid.draw(vd)
grid.text(legend_str_lines[2], x = unit(0.5, "npc"), y = unit(0.04, "npc"),
          just = "center", gp = gpar(fontsize = 20, fontface = "italic"))
dev.off()


cat(
  "✅ Annotated & Venn diagram complete.\n",
  "- Top20_Comp1_annotated.xlsx ->", file.path(output_dir, "Top20_Comp1_annotated.xlsx"), "\n",
  "- SA_up_annotated.xlsx       ->", file.path(output_dir, "Upregulated_in_SA_SA_vs_MOCK_annotated.xlsx"), "\n",
  "- Venn PNG                   ->", venn_out, "\n"
)


################################################################################
################################################################################
################################################################################

# ---- [Export Shared Proteins for Enrichment] ----

# Join top20 and SA-up to retain metadata for shared proteins
shared_annot <- sa_up_annot %>%
  filter(Annotated_Gene %in% shared_proteins) %>%
  distinct(Annotated_Gene, First_ID, Protein_IDs, .keep_all = TRUE)

# Export to Excel for next step
shared_file <- file.path(output_dir, "Shared_Proteins_Top20_vs_SA_vs_MOCK.xlsx")
write.xlsx(shared_annot, file = shared_file, rowNames = FALSE)

cat("\n📤 Exported 11 shared proteins to:\n", shared_file, "\n")

################################################################################
################################################################################
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

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")
output_dir <- "sPLSDA_2D_Plots"

# ---- [Load Shared Proteins] ----
shared_file <- file.path(output_dir, "Shared_Proteins_Top20_vs_SA_vs_MOCK.xlsx")
shared_data <- read.xlsx(shared_file)

# ---- [Prepare Annotation Columns] ----
shared_annot <- shared_data %>%
  mutate(
    First_ID = sub(";.*", "", Protein_IDs),
    Species = case_when(
      grepl("^hsa:|^P\\d|^Q\\d|^O\\d|^A0A", First_ID) ~ "Human",
      grepl("^sav:|^sau:|^SA|^A[0-9]", First_ID, ignore.case = TRUE) ~ "SA",
      TRUE ~ "Unknown"
    )
  )

# ---- [Save Intermediate Annotated File] ----
annot_file <- file.path(output_dir, "Shared_Proteins_annotated_for_enrichment.xlsx")
write.xlsx(shared_annot, annot_file, rowNames = FALSE)

# ---- [Run Enrichment Function] ----
run_enrichment <- function(sig_file, out_dir, label_suffix) {
  if (!file.exists(sig_file)) {
    cat("⚠️ File not found:", sig_file, "\n")
    return()
  }
  
  sig_data <- read.xlsx(sig_file) %>%
    mutate(
      First_ID = sub(";.*", "", Protein_IDs),
      Species = case_when(
        grepl("^hsa:|^P\\d|^Q\\d|^O\\d|^A0A", First_ID) ~ "Human",
        grepl("^sav:|^sau:|^SA|^A[0-9]", First_ID, ignore.case = TRUE) ~ "SA",
        TRUE ~ "Unknown"
      )
    )
  
  human_ids <- gsub("^hsa:", "", sig_data$First_ID[sig_data$Species == "Human"])
  sa_ids <- gsub("^sav:|^sau:", "", sig_data$First_ID[sig_data$Species == "SA"])
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- [Human Enrichment] ----
  valid_keys <- human_ids[human_ids %in% keys(org.Hs.eg.db, keytype = "UNIPROT")]
  
  if (length(valid_keys) > 0) {
    human_entrez <- tryCatch({
      bitr(valid_keys, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) NULL)
    
    if (!is.null(human_entrez) && nrow(human_entrez) > 0) {
      human_genes <- unique(na.omit(human_entrez$ENTREZID))
      
      # GO
      go_enrich <- tryCatch({
        enrichGO(
          gene = human_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = "ALL",
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          readable = TRUE
        )
      }, error = function(e) NULL)
      
      if (!is.null(go_enrich)) {
        write.xlsx(as.data.frame(go_enrich),
                   file = file.path(out_dir, paste0("GO_Human_enrichment_", label_suffix, ".xlsx")),
                   rowNames = FALSE)
      }
      
      # KEGG
      kegg <- tryCatch({
        enrichKEGG(
          gene = human_genes,
          organism = "hsa",
          keyType = "kegg",
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      }, error = function(e) NULL)
      
      if (!is.null(kegg)) {
        kegg_df <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% as.data.frame()
        write.xlsx(kegg_df,
                   file = file.path(out_dir, paste0("KEGG_Human_enrichment_", label_suffix, ".xlsx")),
                   rowNames = FALSE)
      }
      
      # Reactome
      reactome <- tryCatch({
        enrichPathway(
          gene = human_genes,
          organism = "human",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          readable = TRUE
        )
      }, error = function(e) NULL)
      
      if (!is.null(reactome)) {
        reactome_df <- as.data.frame(reactome)
        write.xlsx(reactome_df,
                   file = file.path(out_dir, paste0("Reactome_Human_enrichment_", label_suffix, ".xlsx")),
                   rowNames = FALSE)
      }
    }
  }
  
  # ---- [S. aureus KEGG] ----
  if (length(sa_ids) > 0) {
    sa_kegg <- tryCatch({
      enrichKEGG(
        gene = sa_ids,
        organism = "sau",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
    }, error = function(e) NULL)
    
    if (!is.null(sa_kegg)) {
      sa_df <- as.data.frame(sa_kegg)
      write.xlsx(sa_df,
                 file = file.path(out_dir, paste0("KEGG_SA_enrichment_", label_suffix, ".xlsx")),
                 rowNames = FALSE)
    }
  }
}

# ---- [Run it for Shared Proteins] ----
run_enrichment(
  sig_file = annot_file,
  out_dir = file.path(output_dir, "annotation_shared"),
  label_suffix = "shared_Top20_vs_SA_vs_MOCK"
)

cat("\n✅ Enrichment complete for 11 shared proteins (GO, KEGG, Reactome).\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(readxl)
library(readr)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# ---- [Set Paths] ----
base_dir <- "C:/Users/ga53hil/Desktop/P257_07_new_analysis"
shared_file <- file.path(base_dir, "sPLSDA_2D_Plots", "Shared_Proteins_annotated_for_enrichment.xlsx")
background_file <- file.path(base_dir, "Adjusted_Normalized_P257_07.xlsx")
enrichment_dir <- file.path(base_dir, "sPLSDA_2D_Plots", "annotation_shared")
output_dir <- file.path(base_dir, "sPLSDA_2D_Plots", "Fisher_Test_Shared_Proteins")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Load Data] ----
shared_annot <- read_xlsx(shared_file)
background_data <- read_xlsx(background_file)

# ---- [Prepare Variables] ----
sig_genes <- unique(shared_annot$Annotated_Gene[!is.na(shared_annot$Annotated_Gene)])
background_genes <- unique(background_data$Gene_names)
total_background_size <- length(background_genes)

# ---- [Enrichment Sources] ----
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
suffix <- "shared_Top20_vs_SA_vs_MOCK"

# ---- [Fisher Test Loop] ----
for (source in sources) {
  enrich_file <- file.path(enrichment_dir, paste0(source, "_enrichment_", suffix, ".xlsx"))
  if (!file.exists(enrich_file)) {
    cat("⚠️ Skipping missing file:", enrich_file, "\n")
    next
  }
  
  enrich_data <- read_xlsx(enrich_file)
  
  if (!"geneID" %in% colnames(enrich_data)) {
    cat("⚠️ No geneID column in", enrich_file, "\n")
    next
  }

  # Parse gene lists
  enrich_data <- enrich_data %>%
    mutate(
      gene_list = strsplit(as.character(geneID), "[;/,]"),
      gene_list = lapply(gene_list, trimws)
    )

  # ---- [Fisher's Exact Test] ----
  fisher_results <- data.frame()
  
  for (i in seq_len(nrow(enrich_data))) {
    term <- enrich_data$Description[i]
    term_genes <- intersect(unlist(enrich_data$gene_list[i]), background_genes)
    sig_in_pathway <- intersect(term_genes, sig_genes)
    
    a <- length(sig_in_pathway)
    b <- length(sig_genes) - a
    c <- length(term_genes) - a
    d <- total_background_size - a - b - c
    
    if (any(c(a, b, c, d) < 0)) next
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    p_val <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
    
    if (!is.na(p_val)) {
      gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
      bg_ratio <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
      enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)
      
      fisher_results <- rbind(fisher_results, data.frame(
        Pathway = term,
        Sig_protein_volcano = a,
        Sig_NotIn_Pathway = b,
        Nonsig_In_Pathway = c,
        Nonsig_NotIn_Pathway = d,
        Genes = paste(sig_in_pathway, collapse = ";"),
        p_value = p_val,
        enrichment_factor = enrichment_factor
      ))
    }
  }

  if (nrow(fisher_results) > 0) {
    fisher_results <- fisher_results %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(p_adj)
    
    out_file <- file.path(output_dir, paste0("Fisher_Pathway_", source, "_Shared_Proteins.xlsx"))
    write.xlsx(fisher_results, out_file, rowNames = FALSE)
    cat("✅ Fisher result saved:", out_file, "\n")
  } else {
    cat("⚠️ No valid Fisher results for:", source, "\n")
  }
}

cat("🎉 Fisher’s exact test completed for Shared Proteins.\n")


################################################################################
################################################################################
################################################################################

################################################################################
# Enrichment Heatmap Plotting for Fisher’s Exact Test Results
################################################################################

# ---- [Libraries] ----
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ---- [Set Paths] ----
base_dir <- "C:/Users/ga53hil/Desktop/P257_07_new_analysis"
shared_file <- file.path(base_dir, "sPLSDA_2D_Plots", "Shared_Proteins_annotated_for_enrichment.xlsx")
fisher_dir <- file.path(base_dir, "sPLSDA_2D_Plots", "Fisher_Test_Shared_Proteins")
plot_dir <- file.path(base_dir, "sPLSDA_2D_Plots", "Combined_Human_Enrichment_Heatmap_Sepsis")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---- [Sepsis Keywords] ----
sepsis_keywords <- c("sepsis", "infection", "endotoxin", "bacteria", "pathogen", "defense", "acute", 
                     "positive regulation of innate immune response", "acute inflammation response",
                     "vesicle", "exosome", "positive regulation of exocytosis", "extracellular vesicle",
                     "blood coagulation", "positive regulation of secretion", "blood microparticle",
                     "leukocyte chemotaxis", "lytic vacuole membrane", "secretory granule membrane",
                     "complement activation", "adaptive immune response", 
                     "defense response to Gram-positive bacterium",
                     "positive regulation interleukin-1 production", "macrophage activation",
                     "B cell mediated immunity")

# ---- [Load Shared Annotated Proteins] ----
shared_annot <- read_xlsx(shared_file)
sig_genes <- unique(shared_annot$Annotated_Gene[!is.na(shared_annot$Annotated_Gene)])

# ---- [Sources & File Load] ----
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
combined_all <- list()

for (source in sources) {
  fisher_file <- file.path(fisher_dir, paste0("Fisher_Pathway_", source, "_Shared_Proteins.xlsx"))
  if (!file.exists(fisher_file)) {
    cat("⚠️ Missing:", fisher_file, "\n")
    next
  }
  
  enrich_data <- read.xlsx(fisher_file) %>% as.data.frame()
  if (!"Genes" %in% colnames(enrich_data)) next
  
  enrich_clean <- enrich_data %>%
    filter(!is.na(p_adj), p_adj < 0.05, Sig_protein_volcano > 1, !is.na(Genes), Genes != "") %>%
    mutate(
      Pathway = str_trim(Pathway),
      Source = gsub("_Human", "", source),
      Pathway = paste0(str_wrap(Pathway, 40), " [", Source, "]"),
      log10_p = -log10(p_adj)
    ) %>%
    separate_rows(Genes, sep = "[;,/\\s]+") %>%
    mutate(Gene = str_trim(Genes)) %>%
    filter(Gene != "") %>%
    distinct(Pathway, Gene, log10_p)
  
  combined_all[[source]] <- enrich_clean
}

merged_all <- bind_rows(combined_all)

# ---- [Filter for Sepsis Pathways] ----
sepsis_pathways <- unique(merged_all$Pathway[str_detect(tolower(merged_all$Pathway), paste(tolower(sepsis_keywords), collapse = "|"))])

# Create grid with all shared proteins across all sepsis-related pathways
merged_sepsis <- expand.grid(
  Gene = sig_genes,
  Pathway = sepsis_pathways,
  stringsAsFactors = FALSE
) %>%
  left_join(merged_all, by = c("Gene", "Pathway"))

# ---- [Plotting] ----
if (nrow(merged_sepsis) == 0) {
  cat("⚠️ No sepsis-related enrichment found.\n")
} else {
  # Rank Pathways by average -log10(p.adj)
  pathway_ranking <- merged_sepsis %>%
    group_by(Pathway) %>%
    summarize(avg_log10_p = mean(log10_p, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_log10_p))
  
  merged_sepsis <- merged_sepsis %>%
    mutate(
      Pathway = factor(Pathway, levels = pathway_ranking$Pathway),
      Gene = factor(Gene, levels = sort(sig_genes))
    )
  
  # Split into pages
  chunks <- split(pathway_ranking$Pathway, ceiling(seq_along(pathway_ranking$Pathway) / 20))
  
  for (i in seq_along(chunks)) {
    current_paths <- chunks[[i]]
    chunk_data <- merged_sepsis %>%
      filter(Pathway %in% current_paths) %>%
      mutate(Pathway = factor(Pathway, levels = rev(current_paths)))
    
    plot_title <- "Combined Human Enrichment: shared_sepsis_EVs_proteins"
    plot_file <- file.path(plot_dir, paste0("Combined_Human_Sepsis_Heatmap_Page", i, ".png"))
    
    p <- ggplot(chunk_data, aes(x = Gene, y = Pathway, fill = log10_p)) +
      geom_tile(color = "white", na.rm = FALSE) +
      scale_fill_gradient(
        low = "lightblue", high = "red", na.value = "#4D4D4D", name = "-log10(p.adj)"
      ) +
      labs(
        title = plot_title,
        x = "Shared Significant Protein",
        y = "Sepsis-Enriched Pathway"
      ) +
      theme_minimal(base_size = 18) +  # Increase overall base size
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(plot_file, plot = p,
           width = min(20, max(10, length(unique(chunk_data$Gene)) * 0.3)),
           height = min(50, max(6, 0.4 * length(current_paths))),
           dpi = 300, bg = "white", limitsize = FALSE)
    
    cat("✅ Saved:", plot_file, "\n")
  }
}

cat("\n🎉 Finished: Combined enrichment heatmaps with sepsis filter (shared_sepsis_EVs_proteins).\n")

# ---- [Libraries] ----
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Load Expression Data] ----
data <- read.table("08.12.25_NEW_Adjusted_Normalized_sepsis_P257_09.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# ---- [Load Combined Annotation File] ----
Annotation_Human_SA <- read.table(
  "Annotations(Human+all S. aureus).txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = ""
)
colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))
stopifnot(all(c("Entry", "Gene_Names__primary_") %in% colnames(Annotation_Human_SA)))

# ---- [Clean Gene Names: fallback to Protein_IDs] ----
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "", 
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))

# ---- [Define Expression Groups] ----
group_cols <- list(
  HC = grep("^HC_[1-6]$", colnames(data), value = TRUE),
  BCN = grep("^BCN_[1-6]$", colnames(data), value = TRUE),
  BCP = grep("^BCP_[1-6]$", colnames(data), value = TRUE),
  Bac_infection_group = c(
    grep("^BCN_[1-6]$", colnames(data), value = TRUE),
    grep("^BCP_[1-6]$", colnames(data), value = TRUE)
  )
)

# ---- [Define Comparisons] ----
comparisons <- list(
  "BCN_vs_HC" = c("BCN", "HC"),
  "BCP_vs_HC" = c("BCP", "HC"),
  "BCP_vs_BCN" = c("BCP", "BCN"),
  "Bac_infection_group_vs_HC" = c("Bac_infection_group", "HC")
)

# ---- [Thresholds] ----
pval_threshold <- 0.05
fc_threshold <- 0.3010  # log10(2)

# ---- [Process Each Comparison] ----
for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  label_up <- paste("Higher in", group1)
  label_down <- paste("Higher in", group2)
  
  g1_cols <- group_cols[[group1]]
  g2_cols <- group_cols[[group2]]
  
  dir.create(comp_name, showWarnings = FALSE)
  
  # ---- [Create Data Frame with Expression & IDs] ----
  df <- data.frame(
    Gene = data$Gene_names,
    Protein_IDs = data$Protein_IDs,
    data[, c(g1_cols, g2_cols)],
    stringsAsFactors = FALSE
  )
  
  df$log10FC <- NA
  df$p_value <- NA
  
  for (i in 1:nrow(df)) {
    vals1 <- as.numeric(df[i, g1_cols])
    vals2 <- as.numeric(df[i, g2_cols])
    
    if (all(is.na(vals1)) || all(is.na(vals2))) next
    
    test_result <- tryCatch(t.test(vals1, vals2), error = function(e) NULL)
    if (!is.null(test_result)) {
      df$log10FC[i] <- mean(vals1, na.rm = TRUE) - mean(vals2, na.rm = TRUE)
      df$p_value[i] <- test_result$p.value
    }
  }
  
  df$negLog10P <- -log10(df$p_value)
  
  # ---- [Annotation] ----
  df <- df %>%
    mutate(First_ID = sub(";.*", "", Protein_IDs)) %>%
    left_join(
      Annotation_Human_SA %>%
        select(Entry, Gene_Names__primary_) %>%
        rename(First_ID = Entry, Annotated_Gene = Gene_Names__primary_),
      by = "First_ID"
    )
  
  df$Annotated_Gene[is.na(df$Annotated_Gene) | df$Annotated_Gene == ""] <- 
    df$First_ID[is.na(df$Annotated_Gene) | df$Annotated_Gene == ""]
  
  # ---- [Assign volcano groups] ----
  df$group <- dplyr::case_when(
    df$log10FC > fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_up,
    df$log10FC < -fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_down,
    TRUE ~ "Non-significant"
  )
  
  df$group <- factor(df$group, levels = c(label_up, label_down, "Non-significant"))
  
  # ---- [Export Result Tables] ----
  write.table(df, file = file.path(comp_name, paste0("All_proteins_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df[df$group == label_up, ], 
              file = file.path(comp_name, paste0("Upregulated_in_", group1, "_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df[df$group == label_down, ], 
              file = file.path(comp_name, paste0("Upregulated_in_", group2, "_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # ---- [Top Labels] ----
  top_up <- df %>%
    filter(group == label_up, Annotated_Gene != "") %>%
    arrange(desc(negLog10P), desc(abs(log10FC))) %>%
    head(20)
  
  top_down <- df %>%
    filter(group == label_down, Annotated_Gene != "") %>%
    arrange(desc(negLog10P), desc(abs(log10FC))) %>%
    head(20)
  
  # ---- [Color Map for Plot] ----
  color_map <- setNames(c("red", "blue"), c(label_up, label_down))
  
  # ---- [Volcano Plot] ----
  volcano_plot <- ggplot() +
    geom_point(data = df %>% filter(group == "Non-significant"),
               aes(x = log10FC, y = negLog10P),
               color = "black", alpha = 0.5, size = 2.5, show.legend = FALSE) +
    geom_point(data = df %>% filter(group != "Non-significant"),
               aes(x = log10FC, y = negLog10P, color = group),
               alpha = 0.7, size = 4) +
    scale_color_manual(values = color_map) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black", linewidth = 0.6)
  
  if (nrow(top_up) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = top_up, 
        aes(x = log10FC, y = negLog10P, label = Annotated_Gene),
        color = "red", size = 7, max.overlaps = 100, box.padding = 0.5, segment.size = 0.6, force = 2
      )
  }
  if (nrow(top_down) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = top_down, 
        aes(x = log10FC, y = negLog10P, label = Annotated_Gene),
        color = "blue", size = 7, max.overlaps = 100, box.padding = 0.5, segment.size = 0.6, force = 2
      )
  }
  
  volcano_plot <- volcano_plot +
    labs(
      title = "",
      x = "Log10 Fold Change",
      y = "-Log10 p-value"
    ) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 32),
      axis.title = element_text(size = 28, face = "bold"),
      axis.text = element_text(size = 22),
      legend.position = "bottom",
      legend.box.just = "left", 
      legend.title = element_blank(),
      legend.text = element_text(size = 22),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))
  
  # ---- [Save Plot: Large canvas for big text] ----
  ggsave(
    file.path(comp_name, paste0("Volcano_", comp_name, ".png")),
    plot = volcano_plot, dpi = 300,
    width = 15, height = 12, bg = "white", units = "in"
  )
  
  cat("✅ Volcano plot saved for", comp_name, "\n")
}

cat("🎉 All volcano plots complete with fixed legend and annotation!\n")


################################################################################
################################################################################
################################################################################

# ---- [Step 2: Summary of Output Files] ----
cat("\n📦 Summary of all files generated under P257_09:\n")

comparison_folders <- names(comparisons)  # This ensures consistency with Step 1

for (folder in comparison_folders) {
  comp_path <- file.path(getwd(), folder)
  if (dir.exists(comp_path)) {
    files <- list.files(comp_path, pattern = "\\.(txt|xlsx|png)$", full.names = FALSE)
    cat("\n📁", folder, ":\n")
    if (length(files) == 0) {
      cat("  (No result files found)\n")
    } else {
      print(files)
    }
  } else {
    cat("\n⚠️ Folder missing:", folder, "\n")
  }
}

cat("\n📑 You can now copy these filenames to use in the UpSet plot script.\n")


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
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Use comparisons from Step 1] ----
comparison_folders <- names(comparisons)

# ---- [Enrichment Function] ----
run_enrichment <- function(sig_file, out_dir, label_prefix) {
  if (!file.exists(sig_file)) {
    cat("⚠️ Skipping - file not found:", sig_file, "\n")
    return()
  }
  
  sig_data <- read_tsv(sig_file, show_col_types = FALSE) %>%
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
      go_enrich <- enrichGO(
        gene = human_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
      go_df <- as.data.frame(go_enrich)
      go_df$GeneRatio <- gsub("/", " of ", go_df$GeneRatio)
      go_df$BgRatio <- gsub("/", " of ", go_df$BgRatio)
      write.xlsx(go_df, file = file.path(out_dir, paste0("GO_Human_enrichment_", label_prefix, ".xlsx")), rowNames = FALSE)
      
      # KEGG
      kegg <- enrichKEGG(
        gene = human_genes,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      kegg_df <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% as.data.frame()
      kegg_df$GeneRatio <- gsub("/", " of ", kegg_df$GeneRatio)
      kegg_df$BgRatio <- gsub("/", " of ", kegg_df$BgRatio)
      write.xlsx(kegg_df, file = file.path(out_dir, paste0("KEGG_Human_enrichment_", label_prefix, ".xlsx")), rowNames = FALSE)
      
      # Reactome
      reactome <- enrichPathway(
        gene = human_genes,
        organism = "human",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
      reactome_df <- as.data.frame(reactome)
      reactome_df$GeneRatio <- gsub("/", " of ", reactome_df$GeneRatio)
      reactome_df$BgRatio <- gsub("/", " of ", reactome_df$BgRatio)
      write.xlsx(reactome_df, file = file.path(out_dir, paste0("Reactome_Human_enrichment_", label_prefix, ".xlsx")), rowNames = FALSE)
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
      sa_df$GeneRatio <- gsub("/", " of ", sa_df$GeneRatio)
      sa_df$BgRatio <- gsub("/", " of ", sa_df$BgRatio)
      write.xlsx(sa_df, file = file.path(out_dir, paste0("KEGG_SA_enrichment_", label_prefix, ".xlsx")), rowNames = FALSE)
    }
  }
}

# ---- [Main Loop: Based on Step 1 Comparisons] ----
for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  up_file <- file.path(comp_name, paste0("Upregulated_in_", group1, "_", comp_name, ".txt"))
  down_file <- file.path(comp_name, paste0("Upregulated_in_", group2, "_", comp_name, ".txt"))
  
  up_out <- file.path(comp_name, "annotation_up")
  down_out <- file.path(comp_name, "annotation_down")
  
  cat("\n🔁 Running enrichment for:", comp_name, "\n")
  
  run_enrichment(up_file, up_out, group1)
  run_enrichment(down_file, down_out, group2)
}

cat("\n🎉 Enrichment analysis completed for all comparisons.\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(dplyr)
library(readr)
library(openxlsx)
library(readxl)
library(stringr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Load Data] ----
data <- read.table("08.12.25_NEW_Adjusted_Normalized_sepsis_P257_09.txt",
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

Annotation_Human_SA <- read.table(
  "Annotations(Human+all S. aureus).txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = ""
)
colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))
stopifnot(all(c("Entry", "Gene_Names__primary_") %in% colnames(Annotation_Human_SA)))

# ---- [Base Annotation Matrix] ----
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "",
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))

base_annot <- data %>%
  dplyr::select(Gene_names, Protein_IDs, Majority_protein_IDs,
                Protein_names, Number_of_proteins) %>%
  mutate(First_ID = sub(";.*", "", Protein_IDs)) %>%
  left_join(
    Annotation_Human_SA %>%
      dplyr::select(Entry, Gene_Names__primary_) %>%
      dplyr::rename(First_ID = Entry, Annotated_Gene = Gene_Names__primary_),
    by = "First_ID"
  )

base_annot$Annotated_Gene[is.na(base_annot$Annotated_Gene) | base_annot$Annotated_Gene == ""] <-
  base_annot$First_ID[is.na(base_annot$Annotated_Gene) | base_annot$Annotated_Gene == ""]

# ---- [Use comparisons from Step 1] ----
comparison_names <- names(comparisons)
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")

# ---- [Export Master Annotation for UP and DOWN Proteins] ----
for (comp in comparison_names) {
  cat("\n🧩 Processing:", comp, "\n")
  groups <- comparisons[[comp]]
  up_group <- groups[1]
  down_group <- groups[2]
  
  # Process UP
  up_file <- file.path(comp, paste0("Upregulated_in_", up_group, "_", comp, ".txt"))
  if (file.exists(up_file)) {
    up_data <- read_tsv(up_file, show_col_types = FALSE)
    up_genes <- unique(up_data$Gene)
    
    up_matrix <- base_annot %>%
      mutate(Significant = ifelse(Gene_names %in% up_genes, "+", ""))
    
    out_up <- file.path(comp, "Master_Annotation_up.xlsx")
    write.xlsx(up_matrix, out_up, rowNames = FALSE)
    cat("✅ Exported UP annotation to:", out_up, "\n")
  }
  
  # Process DOWN
  down_file <- file.path(comp, paste0("Upregulated_in_", down_group, "_", comp, ".txt"))
  if (file.exists(down_file)) {
    down_data <- read_tsv(down_file, show_col_types = FALSE)
    down_genes <- unique(down_data$Gene)
    
    down_matrix <- base_annot %>%
      mutate(Significant = ifelse(Gene_names %in% down_genes, "+", ""))
    
    out_down <- file.path(comp, "Master_Annotation_down.xlsx")
    write.xlsx(down_matrix, out_down, rowNames = FALSE)
    cat("✅ Exported DOWN annotation to:", out_down, "\n")
  }
}

# ---- [Fisher Test on Both UP and DOWN] ----
for (comp in comparison_names) {
  cat("\n🔬 Running Fisher's test for:", comp, "\n")
  groups <- comparisons[[comp]]
  up_group <- groups[1]
  down_group <- groups[2]
  
  for (type in c("up", "down")) {
    annot_file <- file.path(comp, paste0("Master_Annotation_", type, ".xlsx"))
    if (!file.exists(annot_file)) {
      cat("⚠️ Missing annotation file for:", type, "\n")
      next
    }
    
    df <- read_xlsx(annot_file)
    df <- df[!is.na(df$Annotated_Gene) & df$Annotated_Gene != "", ]
    background_genes <- unique(df$Annotated_Gene)
    sig_genes <- unique(df$Annotated_Gene[df$Significant == "+"])
    
    fisher_dir <- file.path(comp, paste0("Fisher_", toupper(type)))
    if (!dir.exists(fisher_dir)) dir.create(fisher_dir)
    
    enrich_dir <- file.path(comp, paste0("annotation_", type))
    group_label <- ifelse(type == "up", up_group, down_group)
    
    for (source in sources) {
      enrich_file <- file.path(enrich_dir, paste0(source, "_enrichment_", group_label, ".xlsx"))
      if (!file.exists(enrich_file)) {
        cat("⚠️ Missing enrichment file:", enrich_file, "\n")
        next
      }
      
      enrich_data <- read_xlsx(enrich_file)
      if (!"geneID" %in% colnames(enrich_data)) next
      
      enrich_data <- enrich_data %>%
        mutate(gene_list = strsplit(as.character(geneID), "[,/;]"),
               gene_list = lapply(gene_list, trimws))
      
      fisher_results <- data.frame()
      
      for (i in 1:nrow(enrich_data)) {
        term <- enrich_data$Description[i]
        term_genes <- intersect(unlist(enrich_data$gene_list[i]), background_genes)
        sig_in_pathway <- intersect(term_genes, sig_genes)
        
        a <- length(sig_in_pathway)
        b <- length(sig_genes) - a
        c <- length(term_genes) - a
        d <- length(background_genes) - a - b - c
        
        mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        if (any(mat < 0) || any(!is.finite(mat))) next
        
        test <- tryCatch(fisher.test(mat), error = function(e) NULL)
        if (!is.null(test)) {
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
            p_value = test$p.value,
            enrichment_factor = enrichment_factor
          ))
        }
      }
      
      if (nrow(fisher_results) > 0) {
        fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")
        fisher_results <- fisher_results[order(fisher_results$p_adj), ]
        out_file <- file.path(fisher_dir, paste0("Fisher_Pathway_", source, "_", toupper(type), "_", comp, ".xlsx"))
        write.xlsx(fisher_results, out_file, rowNames = FALSE)
        cat("📁 Saved:", out_file, "\n")
      } else {
        cat("⚠️ No valid Fisher results for", source, type, "in", comp, "\n")
      }
    }
  }
}

cat("\n✅ Fisher's test finished for both UP and DOWN annotations in all comparisons.\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)
library(grid)  # for unit()

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Use comparisons from Step 1] ----
# Ensure 'comparisons' exists in your environment from previous code
comparison_names <- names(comparisons)
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")

# ---- [Sepsis-Related Keywords] ----
sepsis_keywords <- c(
  "sepsis", "infection", "endotoxin", "bacteria",
  "pathogen", "defense", "acute", "positive regulation of innate immune response", 
  "acute inflammation response", "vesicle", "exosome",
  "positive regulation of exocytosis", "extracellular vesicle", "blood coagulation",
  "positive regulation of secrtetion", "blood microparticle", "leukocyte chemotaxis", 
  "lytic vacuole membrane", "secretory granule membrane", "complement activation",
  "positive regulation of adaptive immune response", "defense response to Gram-positive bacterium",
  "positive regulation interleukin-1 production", "macrophage activation", "B cell mediated immunity"
)

exclude_keywords <- c("infections")

# ---- [Generate Enrichment Heatmaps for UP and DOWN] ----
for (comp in comparison_names) {
  for (type in c("UP", "DOWN")) {
    
    plot_folder <- file.path(comp, paste0("Enrichment_Plots_CombinedFiltered_", type))
    if (!dir.exists(plot_folder)) dir.create(plot_folder)
    
    cat("\n📊 Plotting combined filtered heatmaps for:", comp, "(", type, ")\n")
    
    all_enrich <- data.frame()
    
    for (source in sources) {
      fisher_file <- file.path(comp, paste0("Fisher_", type), paste0("Fisher_Pathway_", source, "_", type, "_", comp, ".xlsx"))
      
      if (!file.exists(fisher_file)) {
        cat("⚠️ File not found:", fisher_file, "\n")
        next
      }
      
      enrich_data <- tryCatch({
        read.xlsx(fisher_file) %>% as.data.frame()
      }, error = function(e) {
        cat("❌ Error reading:", fisher_file, "\n")
        return(NULL)
      })
      
      if (is.null(enrich_data) || !"Genes" %in% colnames(enrich_data)) {
        cat("⚠️ Missing required columns in", fisher_file, "\n")
        next
      }
      
      enrich_clean <- enrich_data %>%
        filter(!is.na(p_adj), p_adj < 0.01, Sig_protein_volcano > 1, !is.na(Genes), Genes != "") %>%
        filter(str_detect(tolower(Pathway), paste(tolower(sepsis_keywords), collapse = "|"))) %>%
        filter(!str_detect(tolower(Pathway), paste(tolower(exclude_keywords), collapse = "|"))) %>%
        mutate(
          Pathway = paste0(Pathway, " [", gsub("_Human", "", source), "]"),  # No wrapping, full names
          log10_p = -log10(p_adj),
          Source = source
        ) %>%
        separate_rows(Genes, sep = "[;,/\\s]+") %>%
        mutate(Gene = str_trim(Genes)) %>%
        filter(Gene != "") %>%
        distinct(Pathway, Gene, log10_p, Source)
      
      all_enrich <- bind_rows(all_enrich, enrich_clean)
    }
    
    if (nrow(all_enrich) == 0) {
      cat("⚠️ No combined enrichment data found for", comp, "(", type, ")\n")
      next
    }
    
    # Get all unique genes
    gene_list <- unique(all_enrich$Gene)
    
    # Build full matrix for heatmap
    heatmap_matrix <- all_enrich %>%
      dplyr::select(Pathway, Gene, log10_p) %>%
      tidyr::complete(Pathway, Gene = gene_list)
    
    # Count proteins per pathway for ordering
    pathway_counts <- heatmap_matrix %>%
      filter(!is.na(log10_p)) %>%
      count(Pathway, name = "Gene_Count") %>%
      arrange(desc(Gene_Count))
    
    # Set factor levels to sort by number of proteins
    heatmap_matrix <- heatmap_matrix %>%
      mutate(
        Pathway = factor(Pathway, levels = pathway_counts$Pathway),
        Gene = factor(Gene, levels = sort(unique(Gene)))
      )
    
    # ---- Reasonable chunk size: 15 pathways per plot ----
    ordered_pathways <- levels(heatmap_matrix$Pathway)
    chunks <- split(ordered_pathways, ceiling(seq_along(ordered_pathways) / 30))
    
    for (i in seq_along(chunks)) {
      current_paths <- chunks[[i]]
      chunk_data <- heatmap_matrix %>%
        filter(Pathway %in% current_paths) %>%
        mutate(Pathway = factor(Pathway, levels = rev(current_paths)))  # most enriched on top
      
      # --- AUTOMATIC SIZE ADJUSTMENT ---
      gene_font_size <- 28
      pathway_font_size <- 40
      axis_title_size <- 36
      plot_title_size <- 60
      legend_title_size <- 28
      legend_text_size <- 24
      
      num_genes <- length(unique(chunk_data$Gene))
      # Minimum width 16", maximum 40", scaled by gene count and font size
      width_in <- max(16, min(4 + num_genes * (gene_font_size / 24), 40))
      # Minimum height 9", maximum 25", scaled by number of pathways and font size
      height_in <- max(9, min(25, 1.05 * length(current_paths) * (pathway_font_size / 36)))
      
      heatmap_plot <- ggplot(chunk_data, aes(x = Gene, y = Pathway, fill = log10_p)) +
        geom_tile(color = "white") +
        scale_fill_gradient(
          low = "lightblue",
          high = "red",
          na.value = "#4D4D4D",
          name = expression('-log'[10]*'(p.adj)')
        ) +
        scale_x_discrete(position = "top") +
        labs(
          title = paste("Combined Human Enrichment:", comp, "(Page", i, ")"),
          x = "",
          y = ""
        ) +
        theme_minimal(base_size = 20) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = gene_font_size, face = "bold", color = "black"),
          axis.text.y = element_text(size = pathway_font_size, face = "bold", margin = margin(r = 20), color = "black"),
          axis.title.x = element_text(size = axis_title_size, face = "bold", color = "black", margin = margin(t = 12, b = 8)),
          axis.title.y = element_text(size = axis_title_size, face = "bold", color = "black", margin = margin(r = 16)),
          legend.title = element_text(size = legend_title_size, face = "bold", color = "black"),
          legend.text = element_text(size = legend_text_size, color = "black"),
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(1, "cm"),
          plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5, margin = margin(b = -10), color = "black"),
          plot.title.position = "plot",
          plot.margin = margin(12, 12, 12, 12),
          panel.grid = element_blank(),
          legend.position = "right"
        )
      
      plot_file <- file.path(
        plot_folder,
        paste0("Fisher_Combined_Heatmap_", type, "_", comp, "_Page", i, ".png")
      )
      
      ggsave(
        plot_file,
        plot = heatmap_plot,
        width = width_in,
        height = height_in,
        dpi = 300,
        bg = "white",
        limitsize = FALSE
      )
    }
    
    cat("✅ Combined heatmap(s) saved for", type, "in", comp, "\n")
  }
}

cat("\n🎉 All enrichment heatmaps generated with auto-adjust width/height, NO overlap, large black text, and close title!\n")


################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(UpSetR)
library(dplyr)
library(readr)
library(openxlsx)
library(readxl)
library(tidyr)
library(stringr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Comparisons and Sources] ----
comparisons <- c(
  "BCN_vs_HC",
  "BCP_vs_HC",
  "BCP_vs_BCN"
)
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")

# ---- [Create UpSet Input Lists] ----
up_sets <- list()
down_sets <- list()
all_up_genes <- c()
all_down_genes <- c()

for (comp in comparisons) {
  up_group <- unlist(strsplit(comp, "_vs_"))[2]
  down_group <- unlist(strsplit(comp, "_vs_"))[1]
  
  up_file <- file.path(comp, paste0("Upregulated_in_", up_group, "_", comp, ".txt"))
  down_file <- file.path(comp, paste0("Upregulated_in_", down_group, "_", comp, ".txt"))
  
  if (file.exists(up_file)) {
    up_data <- read_tsv(up_file, show_col_types = FALSE)
    up_sets[[comp]] <- unique(up_data$Gene)
    all_up_genes <- union(all_up_genes, up_sets[[comp]])
  }
  if (file.exists(down_file)) {
    down_data <- read_tsv(down_file, show_col_types = FALSE)
    down_sets[[comp]] <- unique(down_data$Gene)
    all_down_genes <- union(all_down_genes, down_sets[[comp]])
  }
}

# ---- [Build Matrix and Assign Pattern] ----
build_matrix_with_patterns <- function(gene_list, set_list) {
  matrix <- data.frame(Gene = gene_list)
  for (comp in names(set_list)) {
    matrix[[comp]] <- ifelse(matrix$Gene %in% set_list[[comp]], 1, 0)
  }
  matrix$Pattern <- apply(matrix[, -1], 1, paste0, collapse = "")
  matrix
}

up_matrix <- build_matrix_with_patterns(all_up_genes, up_sets)
down_matrix <- build_matrix_with_patterns(all_down_genes, down_sets)

# ---- [Plot UpSet Plots] ----
dir.create("UpSet_Results", showWarnings = FALSE)
png("UpSet_Results/UpSet_Higher_Proteins.png", width = 3000, height = 2000, res = 300)
upset(up_matrix[, names(up_sets)], sets = rev(names(up_sets)), sets.bar.color = "red", main.bar.color = "red",
      mainbar.y.label = "Shared Higher Proteins", sets.x.label = "Proteins per Comparison")
dev.off()

png("UpSet_Results/UpSet_Lower_Proteins.png", width = 3000, height = 2000, res = 300)
upset(down_matrix[, names(down_sets)], sets = rev(names(down_sets)), sets.bar.color = "blue", main.bar.color = "blue",
      mainbar.y.label = "Shared Lower Proteins", sets.x.label = "Proteins per Comparison")
dev.off()

# ---- [Load LFQ and Annotation] ----
data <- read.table("Adjusted_Normalized_remove_Candida_sepsis_P257_09.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "", 
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))
lfq_cols <- grep("^(HC|Bac_culture_).+", colnames(data), value = TRUE)
lfq_data <- data[, c("Gene_names", "Protein_IDs", lfq_cols)]
colnames(lfq_data)[1] <- "Gene"

read_annotation <- function(comp, type) {
  file <- file.path(comp, paste0("Master_Annotation_", type, ".xlsx"))
  if (!file.exists(file)) return(NULL)
  ann <- read.xlsx(file)
  colnames(ann)[colnames(ann) == "Gene_names"] <- "Gene"
  ann[, c("Gene", "Annotated_Gene")]
}

read_pathways <- function(comp, group, type) {
  path_dir <- file.path(comp, paste0("annotation_", type))
  all_paths <- data.frame(Gene = character())
  for (source in sources) {
    f <- file.path(path_dir, paste0(source, "_enrichment_", group, ".xlsx"))
    if (!file.exists(f)) next
    df <- read.xlsx(f)
    if (!"geneID" %in% colnames(df) || !"Description" %in% colnames(df)) next
    df <- df %>% 
      filter(!is.na(geneID)) %>%
      separate_rows(geneID, sep = "[;,/\\s]+") %>%
      mutate(Gene = trimws(geneID)) %>%
      group_by(Gene) %>%
      summarise(Pathways = paste(unique(Description), collapse = "; "), .groups = "drop")
    colnames(df)[2] <- source
    all_paths <- full_join(all_paths, df, by = "Gene")
  }
  all_paths
}

# ---- [Shorten and Clean Sheet Names] ----
shorten_name <- function(name) {
  name <- gsub("BCN", "Bac_neg", name)
  name <- gsub("BCP", "Bac_pos", name)
  name <- gsub("[^A-Za-z0-9_]", "_", name)
  substr(name, 1, 31)
}

# ---- [Write UpSet Pattern Annotation] ----
write_upset_annotation <- function(matrix_df, type, sets) {
  wb <- createWorkbook()
  patterns <- unique(matrix_df$Pattern)
  for (pat in patterns) {
    subset_df <- matrix_df %>% filter(Pattern == pat)
    if (nrow(subset_df) == 0) next
    pattern_bits <- strsplit(pat, "")[[1]]
    active_comps <- names(sets)[pattern_bits == "1"]
    merged <- subset_df
    for (comp in active_comps) {
      group <- ifelse(type == "up", strsplit(comp, "_vs_")[[1]][2], strsplit(comp, "_vs_")[[1]][1])
      ann <- read_annotation(comp, type)
      paths <- read_pathways(comp, group, type)
      if (!is.null(ann)) merged <- left_join(merged, ann, by = "Gene")
      if (!is.null(paths)) merged <- left_join(merged, paths, by = "Gene")
    }
    merged <- left_join(merged, lfq_data, by = "Gene")
    sheet_name <- if (length(active_comps) == 1) {
      shorten_name(active_comps[1])
    } else {
      shorten_name(paste("Shared", paste(active_comps, collapse = "_&_")))
    }
    if (sheet_name %in% names(wb)) {
      sheet_name <- paste0(sheet_name, "_", substr(pat, 1, 5))
    }
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, merged)
  }
  saveWorkbook(wb, file = file.path("UpSet_Results", paste0("UpSet_", type, "_Proteins_Annotated.xlsx")), overwrite = TRUE)
}

# ---- [Run Annotation Export] ----
write_upset_annotation(up_matrix, "up", up_sets)
write_upset_annotation(down_matrix, "down", down_sets)

cat("\n✅ UpSet plots and annotated Excel sheets saved with shortened comparison-based sheet names.\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(UpSetR)
library(dplyr)
library(readr)
library(openxlsx)
library(readxl)
library(tidyr)
library(stringr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Comparisons and Sources] ----
comparisons <- c(
  "BCN_vs_HC",
  "BCP_vs_HC",
  "BCP_vs_BCN"
)
sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")

# ---- [Create UpSet Input Lists] ----
up_sets <- list()
down_sets <- list()
all_genes <- c()
all_sets <- list()

for (comp in comparisons) {
  up_group <- unlist(strsplit(comp, "_vs_"))[2]
  down_group <- unlist(strsplit(comp, "_vs_"))[1]
  
  up_file <- file.path(comp, paste0("Upregulated_in_", up_group, "_", comp, ".txt"))
  down_file <- file.path(comp, paste0("Upregulated_in_", down_group, "_", comp, ".txt"))
  
  if (file.exists(up_file)) {
    up_data <- read_tsv(up_file, show_col_types = FALSE)
    up_sets[[paste0(comp, "_UP")]] <- unique(up_data$Gene)
    all_genes <- union(all_genes, up_sets[[paste0(comp, "_UP")]])
    all_sets[[paste0(comp, "_UP")]] <- unique(up_data$Gene)
  }
  if (file.exists(down_file)) {
    down_data <- read_tsv(down_file, show_col_types = FALSE)
    down_sets[[paste0(comp, "_DOWN")]] <- unique(down_data$Gene)
    all_genes <- union(all_genes, down_sets[[paste0(comp, "_DOWN")]])
    all_sets[[paste0(comp, "_DOWN")]] <- unique(down_data$Gene)
  }
}

# ---- [Build Matrix and Assign Pattern] ----
build_matrix_with_patterns <- function(gene_list, set_list) {
  matrix <- data.frame(Gene = gene_list)
  for (comp in names(set_list)) {
    matrix[[comp]] <- ifelse(matrix$Gene %in% set_list[[comp]], 1, 0)
  }
  matrix$Pattern <- apply(matrix[, -1], 1, paste0, collapse = "")
  matrix
}

combined_matrix <- build_matrix_with_patterns(all_genes, all_sets)

# ---- [Plot Combined UpSet Plot] ----
dir.create("UpSet_Results", showWarnings = FALSE)
png("UpSet_Results/UpSet_All_Significant_Proteins.png", width = 3200, height = 2200, res = 300)
upset(combined_matrix[, names(all_sets)], 
      sets = rev(names(all_sets)), 
      sets.bar.color = "purple", 
      main.bar.color = "purple",
      mainbar.y.label = "Shared Significant Proteins",
      sets.x.label = "Proteins per Comparison")
dev.off()

cat("\n✅ Combined UpSet plot with up and down proteins saved successfully.\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(VennDiagram)
library(grid)
library(dplyr)
library(readr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_09_Sepsis_patients")

# ---- [Comparisons] ----
comparisons <- c(
  "BCN_vs_HC",
  "BCP_vs_HC",
  "BCP_vs_BCN"
)

# ---- [Create Output Directory] ----
venn_dir <- "Venn_Results"
dir.create(venn_dir, showWarnings = FALSE)

# ---- [Initialize Group Proteins List] ----
group_proteins <- list(
  HC = character(),
  BCN = character(),
  BCP = character()
)

# ---- [Collect Unique Significant Proteins by Group] ----
for (comp in comparisons) {
  parts <- unlist(strsplit(comp, "_vs_"))
  group1 <- parts[1]
  group2 <- parts[2]
  
  up_file <- file.path(comp, paste0("Upregulated_in_", group2, "_", comp, ".txt"))
  down_file <- file.path(comp, paste0("Upregulated_in_", group1, "_", comp, ".txt"))
  
  if (file.exists(up_file)) {
    up_data <- read_tsv(up_file, show_col_types = FALSE)
    if (!"Gene" %in% colnames(up_data)) {
      colnames(up_data)[1] <- "Gene"
    }
    group_proteins[[group2]] <- union(group_proteins[[group2]], up_data$Gene)
  } else {
    cat("⚠️ Missing file:", up_file, "\n")
  }
  
  if (file.exists(down_file)) {
    down_data <- read_tsv(down_file, show_col_types = FALSE)
    if (!"Gene" %in% colnames(down_data)) {
      colnames(down_data)[1] <- "Gene"
    }
    group_proteins[[group1]] <- union(group_proteins[[group1]], down_data$Gene)
  } else {
    cat("⚠️ Missing file:", down_file, "\n")
  }
}

# ---- [Extract Group Gene Lists] ----
HC_genes <- unique(group_proteins$HC)
bac_neg_genes <- unique(group_proteins$BCN)
bac_pos_genes <- unique(group_proteins$BCP)

# ---- [Compute Overlaps and Unique Sets] ----
overlap_all     <- Reduce(intersect, list(HC_genes, bac_neg_genes, bac_pos_genes))
unique_HC  <- setdiff(HC_genes, union(bac_neg_genes, bac_pos_genes))
unique_bac_neg  <- setdiff(bac_neg_genes, union(HC_genes, bac_pos_genes))
unique_bac_pos  <- setdiff(bac_pos_genes, union(HC_genes, bac_neg_genes))

overlap_1_2     <- intersect(HC_genes, bac_neg_genes)
overlap_1_3     <- intersect(HC_genes, bac_pos_genes)
overlap_2_3     <- intersect(bac_neg_genes, bac_pos_genes)

# ---- [Write Output Files] ----
write_output <- function(data, name) {
  write.table(data, file.path(venn_dir, paste0(name, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

write_output(overlap_all, "Overlap_All_Groups")
write_output(unique_HC, "Unique_HC")
write_output(unique_bac_neg, "Unique_BCN")
write_output(unique_bac_pos, "Unique_BCP")
write_output(overlap_1_2, "Overlap_HC_BCN")
write_output(overlap_1_3, "Overlap_HC_BCP")
write_output(overlap_2_3, "Overlap_BCN_BCP")

# ---- [Generate Venn Plot Object] ----
venn_plot <- venn.diagram(
  x = list(
    "HC" = HC_genes,
    "BCN" = bac_neg_genes,
    "BCP" = bac_pos_genes
  ),
  category.names = c("HC", "BCN", "BCP"),
  filename = NULL,  # Prevent automatic saving
  output = TRUE,
  fill = c("#F8766D", "#377EB8", "forestgreen"),
  alpha = 0.5,
  cat.col = c("#F8766D", "#377EB8", "forestgreen"),
  cat.cex = 1.8,
  cex = 2.5,
  cat.pos = c(-15, 10, 180),
  cat.dist = c(0.05, 0.05, 0.04),
  main = "Venn Diagram: Significant Proteins",
  main.cex = 2,
  euler.d = FALSE,
  scaled = FALSE
)

# ---- [Display in RStudio] ----
grid.newpage()
grid.draw(venn_plot)

# ---- [Save to File as PNG] ----
png(file.path(venn_dir, "Venn_Diagram_Significant_Proteins.png"), width = 2000, height = 2000, res = 300)
grid.draw(venn_plot)
dev.off()

cat("\n✅ Venn diagram displayed and saved in 'Venn_Results' folder.\n")

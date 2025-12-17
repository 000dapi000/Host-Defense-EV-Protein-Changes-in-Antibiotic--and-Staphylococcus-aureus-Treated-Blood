# ---- [Libraries] ----
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)

# ---- [Set Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")

# ---- [Load Data] ----
data <- read.table("Adjusted_Normalized_P257_07.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# ---- [Load Annotation File] ----
Annotation_Human_SA <- read.table(
  "Annotations(Human+all S. aureus).txt",
  header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = ""
)
colnames(Annotation_Human_SA) <- gsub("[[:punct:]]", "_", colnames(Annotation_Human_SA))
stopifnot(all(c("Entry", "Gene_Names__primary_") %in% colnames(Annotation_Human_SA)))

# ---- [Gene Cleanup] ----
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "", 
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))

# ---- [Thresholds] ----
pval_threshold <- 0.05
fc_threshold <- 0.3010  # log10(2)

# ---- [Manual Comparisons] ----
comparisons <- list(
  list(name = "SA_vs_MOCK", group1_name = "SA", group2_name = "MOCK",
       group1_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6"),
       group2_cols = c("MOCK_1", "MOCK_2", "MOCK_3", "MOCK_4", "MOCK_5", "MOCK_6")),
  list(name = "Pip_Tazo_vs_SA", group1_name = "Pip_Tazo", group2_name = "SA",
       group1_cols = c("Pip_Tazo_1", "Pip_Tazo_2", "Pip_Tazo_3", "Pip_Tazo_4", "Pip_Tazo_5", "Pip_Tazo_6"),
       group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6")),
  list(name = "Vancomycin_vs_SA", group1_name = "Vancomycin", group2_name = "SA",
       group1_cols = c("Vancomycin_1", "Vancomycin_2", "Vancomycin_3", "Vancomycin_4", "Vancomycin_5", "Vancomycin_6"),
       group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6")),
  list(name = "Moxifloxacin_vs_SA", group1_name = "Moxifloxacin", group2_name = "SA",
       group1_cols = c("Moxifloxacin_1", "Moxifloxacin_2", "Moxifloxacin_3", "Moxifloxacin_4", "Moxifloxacin_5", "Moxifloxacin_6"),
       group2_cols = c("SA_1", "SA_2", "SA_3", "SA_4", "SA_5", "SA_6"))
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
  
  dir.create(comp_name, showWarnings = FALSE)
  
  df <- data.frame(
    Gene = data$Gene_names,
    Protein_IDs = data$Protein_IDs,
    data[, c(g1_cols, g2_cols)],
    stringsAsFactors = FALSE
  )
  
  df$log10FC <- NA
  df$p_value <- NA
  
  for (i in 1:nrow(df)) {
    vals1 <- as.numeric(unlist(df[i, g1_cols]))
    vals2 <- as.numeric(unlist(df[i, g2_cols]))
    
    if (all(is.na(vals1)) || all(is.na(vals2))) next
    
    test_result <- tryCatch(t.test(vals1, vals2, paired = TRUE), error = function(e) NULL)
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
        dplyr::select(Entry, Gene_Names__primary_) %>%
        dplyr::rename(First_ID = Entry, Annotated_Gene = Gene_Names__primary_),
      by = "First_ID"
    )
  
  df$Annotated_Gene[is.na(df$Annotated_Gene) | df$Annotated_Gene == ""] <- 
    df$First_ID[is.na(df$Annotated_Gene) | df$Annotated_Gene == ""]
  
  # ---- [Volcano Grouping] ----
  df$group <- case_when(
    df$log10FC > fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_up,
    df$log10FC < -fc_threshold & df$negLog10P > -log10(pval_threshold) ~ label_down,
    TRUE ~ "Non-significant"
  )
  df$group <- factor(df$group, levels = c(label_up, label_down, "Non-significant"))
  
  # ---- [Top Genes] ----
  top_up <- df %>%
    filter(group == label_up & !is.na(Annotated_Gene)) %>%
    arrange(desc(negLog10P), desc(abs(log10FC))) %>%
    head(15)
  
  top_down <- df %>%
    filter(group == label_down & !is.na(Annotated_Gene)) %>%
    arrange(desc(negLog10P), desc(abs(log10FC))) %>%
    head(15)
  
  # ---- [Color Map] ----
  color_map <- setNames(c("red", "blue", "black"),
                        c(label_up, label_down, "Non-significant"))
  
  # ---- [Save Up/Downregulated Files] ----
  write.table(df, file = file.path(comp_name, paste0("All_proteins_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df[df$group == label_up, ],
              file = file.path(comp_name, paste0("Upregulated_in_", group1, "_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df[df$group == label_down, ],
              file = file.path(comp_name, paste0("Upregulated_in_", group2, "_", comp_name, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # ---- [Plot Volcano] ----
  volcano_plot <- ggplot(df, aes(x = log10FC, y = negLog10P)) +
    geom_point(aes(color = group), alpha = 0.6, size = 5) +
    scale_color_manual(values = color_map) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
    geom_text_repel(data = top_up, aes(label = Annotated_Gene), color = "red", size = 10) +
    geom_text_repel(data = top_down, aes(label = Annotated_Gene), color = "blue", size = 10) +
    labs(
      title = "Volcano Plot",
      x = "Log10 Fold Change",
      y = "-Log10 p-value"
    ) +
    theme_minimal(base_size = 25) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 25),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 6)))
  
  
  ggsave(file.path(comp_name, paste0("Volcano_", comp_name, ".png")),
         plot = volcano_plot, dpi = 300, width = 12, height = 12, bg = "white")
  
  print(volcano_plot)
  cat("✅ Volcano plot saved for", comp_name, "\n")
}

cat("🎉 All volcano plots and result tables completed for P257_07!\n")


################################################################################
################################################################################
################################################################################

# ---- [Improved Summary of Comparison Files for P257_07] ----
cat("\n📦 Summary of all files generated under P257_07:\n")

comparison_folders <- sapply(comparisons, function(comp) comp$name)

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
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")

# ---- [Comparisons List] ----
comparisons <- list(
  list(name = "SA_vs_MOCK", group1 = "SA", group2 = "MOCK"),
  list(name = "Pip_Tazo_vs_SA", group1 = "Pip_Tazo", group2 = "SA"),
  list(name = "Vancomycin_vs_SA", group1 = "Vancomycin", group2 = "SA"),
  list(name = "Moxifloxacin_vs_SA", group1 = "Moxifloxacin", group2 = "SA")
)

# ---- [Enrichment Function] ----
run_enrichment <- function(sig_file, out_dir, label_suffix) {
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
        write.xlsx(as.data.frame(go_enrich), file = file.path(out_dir, paste0("GO_Human_enrichment_", label_suffix, ".xlsx")), rowNames = FALSE)
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
        write.xlsx(kegg_df, file = file.path(out_dir, paste0("KEGG_Human_enrichment_", label_suffix, ".xlsx")), rowNames = FALSE)
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
        write.xlsx(reactome_df, file = file.path(out_dir, paste0("Reactome_Human_enrichment_", label_suffix, ".xlsx")), rowNames = FALSE)
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
      write.xlsx(sa_df, file = file.path(out_dir, paste0("KEGG_SA_enrichment_", label_suffix, ".xlsx")), rowNames = FALSE)
    }
  }
}

# ---- [Main Loop for All Comparisons] ----
for (comp in comparisons) {
  comp_name <- comp$name
  group1 <- comp$group1
  group2 <- comp$group2
  
  # Corrected File paths (matching updated Step 1)
  up_file <- file.path(comp_name, paste0("Upregulated_in_", group1, "_", comp_name, ".txt"))
  down_file <- file.path(comp_name, paste0("Upregulated_in_", group2, "_", comp_name, ".txt"))
  sig_file <- file.path(comp_name, paste0("significant_proteins_", comp_name, ".txt"))
  
  # Output dirs
  up_out <- file.path(comp_name, "annotation_up")
  down_out <- file.path(comp_name, "annotation_down")
  sig_out <- file.path(comp_name, "annotation_significant")
  
  # Labels
  up_label <- paste0("up_", comp_name)
  down_label <- paste0("down_", comp_name)
  sig_label <- paste0("significant_", comp_name)
  
  cat("\n🔁 Processing:", comp_name, "\n")
  
  run_enrichment(up_file, up_out, up_label)
  run_enrichment(down_file, down_out, down_label)
  run_enrichment(sig_file, sig_out, sig_label)
}

cat("\n🎉 All enrichment analyses completed successfully and saved.\n")

################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(dplyr)
library(readr)
library(readxl)
library(openxlsx)
library(stringr)
library(tidyr)

# ---- [Working Directory] ----
setwd("C:/Users/ga53hil/Desktop/P257_07_new_analysis")

# ---- [Load Base Data] ----
data <- read.table("Adjusted_Normalized_P257_07.txt",
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

required_cols <- c("Entry", "Gene_Names__primary_", "Organism",
                   "Gene_Ontology__GO_", "KEGG", "Reactome")
missing_cols <- setdiff(required_cols, colnames(Annotation_Human_SA))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# ---- [Base Annotation Table with SA Pathways] ----
data$Gene_names <- ifelse(is.na(data$Gene_names) | trimws(data$Gene_names) == "",
                          data$Protein_IDs, data$Gene_names)
data$Gene_names <- make.unique(as.character(data$Gene_names))
data <- data %>% mutate(First_ID = sub(";.*", "", Protein_IDs))

base_annot <- data %>%
  dplyr::select(Gene_names, Protein_IDs, First_ID) %>%
  left_join(
    Annotation_Human_SA %>%
      dplyr::select(Entry, Gene_Names__primary_) %>%
      dplyr::rename(First_ID = Entry, Annotated_Gene = Gene_Names__primary_),
    by = "First_ID"
  )

base_annot$Annotated_Gene[is.na(base_annot$Annotated_Gene) | base_annot$Annotated_Gene == ""] <-
  base_annot$First_ID[is.na(base_annot$Annotated_Gene) | base_annot$Annotated_Gene == ""]

SA_annotation <- Annotation_Human_SA %>%
  filter(!grepl("Homo sapiens", Organism, ignore.case = TRUE)) %>%
  dplyr::select(
    Entry,
    SA_GO = Gene_Ontology__GO_,
    SA_KEGG = KEGG,
    SA_Reactome = Reactome
  )

base_annot <- base_annot %>%
  left_join(SA_annotation, by = c("First_ID" = "Entry"))

write.xlsx(base_annot, "Base_Annotation_with_SA_Pathways.xlsx", rowNames = FALSE)

# ---- [Comparisons and Sources] ----
comparisons <- list(
  SA_vs_MOCK = list(group1 = "SA", group2 = "MOCK"),
  Pip_Tazo_vs_SA = list(group1 = "Pip_Tazo", group2 = "SA"),
  Vancomycin_vs_SA = list(group1 = "Vancomycin", group2 = "SA"),
  Moxifloxacin_vs_SA = list(group1 = "Moxifloxacin", group2 = "SA")
)

sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
sa_sources <- c("SA_GO", "SA_KEGG", "SA_Reactome")

# ---- [Export Master Annotations for UP and DOWN] ----
for (comp in names(comparisons)) {
  cat("\n🧩 Exporting Master Annotations:", comp, "\n")
  
  group1 <- comparisons[[comp]]$group1
  group2 <- comparisons[[comp]]$group2
  
  # UP (higher in group1)
  up_file <- file.path(comp, paste0("Upregulated_in_", group1, "_", comp, ".txt"))
  if (file.exists(up_file)) {
    up_data <- read_tsv(up_file, show_col_types = FALSE)
    up_genes <- unique(up_data$Gene)
    
    up_matrix <- base_annot %>%
      mutate(Significant = ifelse(Gene_names %in% up_genes, "+", ""))
    
    write.xlsx(up_matrix, file.path(comp, "Master_Annotation_up.xlsx"), rowNames = FALSE)
  }
  
  # DOWN (higher in group2)
  down_file <- file.path(comp, paste0("Upregulated_in_", group2, "_", comp, ".txt"))
  if (file.exists(down_file)) {
    down_data <- read_tsv(down_file, show_col_types = FALSE)
    down_genes <- unique(down_data$Gene)
    
    down_matrix <- base_annot %>%
      mutate(Significant = ifelse(Gene_names %in% down_genes, "+", ""))
    
    write.xlsx(down_matrix, file.path(comp, "Master_Annotation_down.xlsx"), rowNames = FALSE)
  }
}

# ---- [Fisher's Exact Test for Enriched Pathways] ----
for (comp in names(comparisons)) {
  cat("\n🔬 Fisher’s Test for:", comp, "\n")
  
  for (type in c("up", "down")) {
    annot_file <- file.path(comp, paste0("Master_Annotation_", type, ".xlsx"))
    if (!file.exists(annot_file)) {
      cat("⚠️ Missing annotation file for", type, "in", comp, "\n")
      next
    }
    
    df <- read_xlsx(annot_file)
    df <- df[!is.na(df$Annotated_Gene) & df$Annotated_Gene != "", ]
    background_genes <- unique(df$Annotated_Gene)
    sig_genes <- unique(df$Annotated_Gene[df$Significant == "+"])
    
    fisher_dir <- file.path(comp, paste0("Fisher_", toupper(type)))
    if (!dir.exists(fisher_dir)) dir.create(fisher_dir)
    
    enrich_dir <- file.path(comp, paste0("annotation_", type))
    label_suffix <- paste0(type, "_", comp)
    
    # ---- [Human Enrichment Sources] ----
    for (source in sources) {
      enrich_file <- file.path(enrich_dir, paste0(source, "_enrichment_", label_suffix, ".xlsx"))
      cat("🔎 Looking for:", enrich_file, "\n")
      
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
      }
    }
    
    # ---- [Fisher for SA-specific annotation columns] ----
    for (sa_source in sa_sources) {
      if (!sa_source %in% colnames(df)) next
      
      sa_results <- df %>%
        filter(!is.na(.data[[sa_source]]), .data[[sa_source]] != "") %>%
        mutate(Pathway = strsplit(as.character(.data[[sa_source]]), "[;,|]")) %>%
        unnest(Pathway) %>%
        mutate(Pathway = str_trim(Pathway)) %>%
        filter(Pathway != "") %>%
        group_by(Pathway) %>%
        summarize(
          Genes = paste(unique(Annotated_Gene[Significant == "+"]), collapse = ";"),
          Sig_protein_volcano = sum(Significant == "+" & Annotated_Gene != ""),
          Nonsig_In_Pathway = sum(Significant != "+" & Annotated_Gene != ""),
          .groups = "drop"
        ) %>%
        rowwise() %>%
        mutate(
          Sig_NotIn_Pathway = length(sig_genes) - Sig_protein_volcano,
          Nonsig_NotIn_Pathway = length(background_genes) - Sig_protein_volcano - Nonsig_In_Pathway - Sig_NotIn_Pathway,
          p_value = tryCatch({
            fisher.test(matrix(c(Sig_protein_volcano, Sig_NotIn_Pathway,
                                 Nonsig_In_Pathway, Nonsig_NotIn_Pathway), nrow = 2))$p.value
          }, error = function(e) NA_real_)
        ) %>%
        ungroup() %>%
        filter(!is.na(p_value)) %>%
        mutate(
          enrichment_factor = (Sig_protein_volcano / length(sig_genes)) /
            ((Sig_protein_volcano + Nonsig_In_Pathway) / length(background_genes)),
          p_adj = p.adjust(p_value, method = "BH")
        ) %>%
        arrange(p_adj)
      
      if (nrow(sa_results) > 0) {
        out_file <- file.path(fisher_dir, paste0("Fisher_Pathway_", sa_source, "_", toupper(type), "_", comp, ".xlsx"))
        write.xlsx(sa_results, out_file, rowNames = FALSE)
        cat("✅ SA Fisher saved:", out_file, "\n")
      }
    }
  }
}

cat("\n✅ Fisher’s exact test (Human + SA) completed for all comparisons.\n")


################################################################################
################################################################################
################################################################################

# ---- [Libraries] ----
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)

# ---- [Sepsis Keywords] ----
sepsis_keywords <- c(
  "sepsis", "infection", "endotoxin", "bacteria",
  "pathogen", "defense", "acute", 
  "positive regulation of innate immune response", 
  "acute inflammation response", "vesicle", 
  "exosome", "positive regulation of exocytosis", 
  "extracellular vesicle", "blood coagulation", 
  "positive regulation of secrtetion", "blood microparticle", 
  "leukocyte chemotaxis", "lytic vacuole membrane", 
  "secretory granule membrane", "complement activation",
  "adaptive immune response", 
  "defense response to Gram-positive bacterium",
  "positive regulation interleukin-1 production", 
  "macrophage activation", "B cell mediated immunity"
)

# ---- [Exclude Non-Sepsis Keywords] ----
exclude_keywords <- c(
  "development", "morphogenesis", "neuron", "axon", "synapse",
  "brain", "behavior", "eye", "ear", "photoreceptor", "limb", "skeletal"
)

# ---- [Merged Human Enrichment Heatmap: SA_vs_MOCK - UP only] ----
cat("\n🧬 Merging GO_Human, KEGG_Human, Reactome_Human into one plot for SA_vs_MOCK (UP only)...\n")

comp <- "SA_vs_MOCK"
type <- "UP"
human_sources <- c("GO_Human", "KEGG_Human", "Reactome_Human")
combined_data <- list()

for (source in human_sources) {
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
    filter(!is.na(p_adj), p_adj < 0.05, Sig_protein_volcano > 2, !is.na(Genes), Genes != "") %>%
    filter(str_detect(tolower(Pathway), paste(tolower(sepsis_keywords), collapse = "|"))) %>%
    filter(!str_detect(tolower(Pathway), paste(tolower(exclude_keywords), collapse = "|"))) %>%
    mutate(
      Pathway = str_trim(Pathway),
      Source = gsub("_Human", "", source),
      Pathway = paste0(str_wrap(Pathway, width = 40), " [", Source, "]"),
      log10_p = -log10(p_adj)
    ) %>%
    separate_rows(Genes, sep = "[;,/\\s]+") %>%
    mutate(Gene = str_trim(Genes)) %>%
    filter(Gene != "") %>%
    distinct(Pathway, Gene, log10_p)
  
  combined_data[[source]] <- enrich_clean
}

combined_all <- bind_rows(combined_data)

if (nrow(combined_all) == 0) {
  cat("⚠️ No enrichment terms found for merging.\n")
} else {
  # Build full matrix
  gene_list <- unique(combined_all$Gene)
  full_matrix <- combined_all %>%
    dplyr::select(Pathway, Gene, log10_p) %>%
    tidyr::complete(Pathway, Gene = gene_list)
  
  # Count and sort pathways by gene count (top = most enriched)
  pathway_counts <- full_matrix %>%
    filter(!is.na(log10_p)) %>%
    count(Pathway, name = "Gene_Count") %>%
    arrange(desc(Gene_Count))
  
  full_matrix <- full_matrix %>%
    mutate(
      Pathway = factor(Pathway, levels = pathway_counts$Pathway),
      Gene = factor(Gene, levels = sort(unique(Gene)))
    )
  
  ordered_pathways <- levels(full_matrix$Pathway)
  chunks <- split(ordered_pathways, ceiling(seq_along(ordered_pathways) / 20))
  
  plot_folder <- file.path(comp, "Combined_Enrichment_Plots_UP")
  if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
  
  for (i in seq_along(chunks)) {
    current_paths <- chunks[[i]]
    chunk_data <- full_matrix %>%
      filter(Pathway %in% current_paths) %>%
      mutate(Pathway = factor(Pathway, levels = rev(current_paths)))  # most enriched on top
    
    plot_title <- paste0("Combined Human Enrichment (GO + KEGG + Reactome): UP - ", comp,
                         " (Sepsis-relevant only, Page ", i, ")")
    
    plot_file <- file.path(
      plot_folder,
      paste0("Combined_Heatmap_", comp, "_UP_Sepsis_Page", i, ".png")
    )
    
    heatmap_plot <- ggplot(chunk_data, aes(x = Gene, y = Pathway, fill = log10_p)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "lightblue", high = "red", na.value = "#4D4D4D", name = "-log10(p.adj)") +
      scale_x_discrete(position = "top") +
      labs(
        title = plot_title,
        x = "Shared Significant Protein",
        y = "Sepsis-Enriched Pathway"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, lineheight = 1.1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      )
    
    ggsave(plot_file,
           plot = heatmap_plot,
           width = min(30, max(12, length(unique(chunk_data$Gene)) * 0.3)),
           height = min(50, max(6, 0.5 * length(current_paths))),
           dpi = 300, bg = "white", limitsize = FALSE)
  }
  
  cat("✅ Combined Human heatmap generated for SA_vs_MOCK - UP - Sepsis\n")
}

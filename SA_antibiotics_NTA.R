# ===============================
# Libraries
# ===============================
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggpubr)

# ===============================
# Output directory
# ===============================
out_dir <- "NTA_Antibiotics_Output"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ===============================
# Read data
# ===============================
data <- read.table(
  file = "C:/Users/ga53hil/Desktop/NTA_antibitoics_SA.txt",
  header = TRUE,
  sep = "\t",
  fileEncoding = "UTF-16LE",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ===============================
# Rename column safely
# ===============================
data <- data %>%
  rename(
    CMO_particles = `CMO_Particles_in _1_ml_plasma`
  )

# ===============================
# DEFINE VALID GROUP NAMES (EXACT)
# ===============================
group_levels <- c(
  "MOCK",
  "SA",
  "Pip-Tazo 25 μg/ml",
  "Pip-Tazo 50 μg/ml",
  "Pip-Tazo 100 μg/ml",
  "Vancomycin 6.25 μg/ml",
  "Vancomycin 12.5 μg/ml",
  "Vancomycin 25 μg/ml",
  "Moxifloxacin 0.18125 μg/ml",
  "Moxifloxacin 0.725 μg/ml",
  "Moxifloxacin 2.9 μg/ml"
)

# ===============================
# Assign GROUP using grepl (SAFE)
# ===============================
data$Group <- NA_character_

for (g in group_levels) {
  data$Group[grepl(paste0("^", g, "_"), data$Samples)] <- g
}

# ===============================
# STOP if any group is unmatched
# ===============================
if (any(is.na(data$Group))) {
  stop(
    "ERROR: Some Samples do not match any predefined group:\n",
    paste(unique(data$Samples[is.na(data$Group)]), collapse = "\n")
  )
}

# ===============================
# Extract PAIR (SAFE)
# ===============================
data$Pair <- sub(".*_", "", data$Samples)

# ===============================
# Convert Group to factor (ORDER ONLY FOR DISPLAY)
# ===============================
data$Group <- factor(data$Group, levels = group_levels)

# ===============================
# Colors
# ===============================
group_colors <- c(
  "MOCK"                       = "#4A4A4A",
  "SA"                         = "#7B3294",
  "Pip-Tazo 25 μg/ml"          = "#9ADBE8",
  "Pip-Tazo 50 μg/ml"          = "#4DBBD5",
  "Pip-Tazo 100 μg/ml"         = "#0085AD",
  "Vancomycin 6.25 μg/ml"      = "#F4A582",
  "Vancomycin 12.5 μg/ml"      = "#E64B35",
  "Vancomycin 25 μg/ml"        = "#B2182B",
  "Moxifloxacin 0.18125 μg/ml" = "#A6DBA0",
  "Moxifloxacin 0.725 μg/ml"   = "#00A087",
  "Moxifloxacin 2.9 μg/ml"     = "#006837"
)

# ===============================
# Publication theme (100% BLACK)
# ===============================
theme_pub_black <- theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    axis.line.x  = element_line(size = 2, color = "black"),
    axis.line.y  = element_line(size = 2, color = "black"),
    axis.ticks.x = element_line(size = 2, color = "black"),
    axis.ticks.y = element_line(size = 2, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    plot.title = element_text(size = 26, hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

# ===============================
# Force Y-axis to start at 0
# ===============================
y_scale_zero <- scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0, 0.05))
)

# ===============================
# Paired Wilcoxon – ALL 1 vs 1
# ===============================
pairwise_wilcoxon_paired <- function(df, value_col) {
  
  combs <- combn(levels(df$Group), 2, simplify = FALSE)
  
  map_dfr(combs, function(g) {
    
    sub <- df %>%
      filter(Group %in% g) %>%
      select(Pair, Group, value = !!sym(value_col)) %>%
      pivot_wider(names_from = Group, values_from = value) %>%
      drop_na()
    
    if (nrow(sub) < 2) return(NULL)
    
    w <- wilcox.test(
      sub[[g[1]]],
      sub[[g[2]]],
      paired = TRUE,
      exact = FALSE
    )
    
    tibble(
      Group1  = g[1],
      Group2  = g[2],
      n_pairs = nrow(sub),
      p_value = w$p.value
    )
  }) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

# ===============================
# Run statistics
# ===============================
wilcox_size <- pairwise_wilcoxon_paired(data, "Particle_size_nm")
wilcox_conc <- pairwise_wilcoxon_paired(data, "CMO_particles")

# ===============================
# Export TSV files
# ===============================
write.table(
  wilcox_size,
  file.path(out_dir, "Wilcoxon_Paired_Particle_Size_All_Pairwise.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  wilcox_conc,
  file.path(out_dir, "Wilcoxon_Paired_CMO_Concentration_All_Pairwise.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ===============================
# Plot: Particle Size
# ===============================
p1 <- ggplot(data, aes(Group, Particle_size_nm, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, size = 1.4, color = "black") +
  geom_jitter(width = 0.12, size = 3.5, shape = 21,
              color = "black", stroke = 1.3) +
  scale_fill_manual(values = group_colors) +
  y_scale_zero +
  labs(
    x = "",
    y = "Particle Size (nm)",
    title = "Particle Size",
    subtitle = "Wilcoxon matched-pairs signed-rank test (all pairwise comparisons)"
  ) +
  theme_pub_black +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# ===============================
# Plot: CMO concentration
# ===============================
p2 <- ggplot(data, aes(Group, CMO_particles, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, size = 1.4, color = "black") +
  geom_jitter(width = 0.12, size = 3.5, shape = 21,
              color = "black", stroke = 1.3) +
  scale_fill_manual(values = group_colors) +
  y_scale_zero +
  labs(
    x = "",
    y = "CMO+ Particles per mL Plasma",
    title = "CMO+ Particle Concentration",
    subtitle = "Wilcoxon matched-pairs signed-rank test (all pairwise comparisons)"
  ) +
  theme_pub_black +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
# ===============================
# Save plots
# ===============================
ggsave(
  file.path(out_dir, "Particle_Size_bar_dots.png"),
  p1, width = 10, height = 7, dpi = 300, bg = "white"
)

ggsave(
  file.path(out_dir, "CMO_Concentration_bar_dots.png"),
  p2, width = 10, height = 7, dpi = 300, bg = "white"
)

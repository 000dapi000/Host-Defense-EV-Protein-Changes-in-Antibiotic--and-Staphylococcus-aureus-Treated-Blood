# ===============================
# Libraries
# ===============================
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggpubr)
library(scales)

# ===============================
# Output directory
# ===============================
out_dir <- "EV_flow_Output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===============================
# Read data
# ===============================
data <- read.table(
  file = "C:/Users/ga53hil/Desktop/15.12.25 EV BacB flow result.txt",
  header = TRUE,
  sep = "\t",
  fileEncoding = "UTF-16LE",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ===============================
# Define group names (display order)
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
# Assign GROUP (case-insensitive)
# ===============================
data$Group <- NA_character_
for (g in group_levels) {
  data$Group[grepl(paste0("^", g, "_"), data$Samples, ignore.case = TRUE)] <- g
}
if (any(is.na(data$Group))) {
  stop(
    "ERROR: Unmatched samples:\n",
    paste(unique(data$Samples[is.na(data$Group)]), collapse = "\n")
  )
}

# ===============================
# Extract pair ID
# ===============================
data$Pair <- sub(".*_", "", data$Samples)

# ===============================
# Factor order
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
# Publication theme (black axes)
# ===============================
theme_force_black_axes <- theme(
  axis.line  = element_line(color = "black", linewidth = 2.2),
  axis.ticks = element_line(color = "black", linewidth = 2.2),
  axis.text  = element_text(color = "black"),
  axis.title = element_text(color = "black")
)

theme_pub_black <- theme_classic(base_size = 22) +
  theme_force_black_axes +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0.25, "cm")
  )

y_scale_zero <- scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0, 0.05))
)

# ===============================
# Paired Wilcoxon (export only)
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
      Group1 = g[1],
      Group2 = g[2],
      n_pairs = nrow(sub),
      p_value = w$p.value
    )
  }) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

wilcox_bmv <- pairwise_wilcoxon_paired(
  data,
  "Percentage of CMD+ PanEV- Anti-SA+ of bMVs-BacB complex"
)

wilcox_ev <- pairwise_wilcoxon_paired(
  data,
  "Percentage of CMD+ PanEV+ Anti-SA+ of EV-Miltenyi beads complex"
)

write.table(
  wilcox_bmv,
  file.path(out_dir, "Wilcoxon_CMD_PanEVneg_bMVs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  wilcox_ev,
  file.path(out_dir, "Wilcoxon_CMD_PanEVpos_EV.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ===============================
# Manual significance tables (USED IN PLOTS)
#   IMPORTANT: stat_pvalue_manual needs y.position
# ===============================
sig_bmv <- tibble(
  group1 = "MOCK",
  group2 = "SA",
  label  = "*",
  y.position = 9
)

sig_miltenyi <- tibble(
  group1 = c("MOCK", "SA", "SA"),
  group2 = c("SA", "Pip-Tazo 25 μg/ml", "Vancomycin 25 μg/ml"),
  label  = "*",
  y.position = max(
    data$`Percentage of CMD+ PanEV+ Anti-SA+ of EV-Miltenyi beads complex`,
    na.rm = TRUE
  ) * c(1.10, 1.18, 1.26)
)

# ===============================
# Piecewise linear transformation
# ===============================
piecewise_linear_trans <- trans_new(
  name = "piecewise_linear",
  transform = function(y) {
    ifelse(y <= 0.5, y * 4, 2 + log10(y + 1) * 4)
  },
  inverse = function(y) {
    ifelse(y <= 2, y / 4, 10^((y - 2) / 4) - 1)
  }
)

# ===============================
# Plot: bMVs–BacB Complex (FIXED significance)
# ===============================
p1 <- ggplot(
  data,
  aes(
    x = Group,
    y = `Percentage of CMD+ PanEV- Anti-SA+ of bMVs-BacB complex`,
    fill = Group
  )
) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 1.5, color = "black") +
  geom_jitter(width = 0.12, size = 3.6, shape = 21,
              color = "black", stroke = 1.3) +
  stat_pvalue_manual(
    sig_bmv,
    tip.length = 0.02,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    trans = piecewise_linear_trans,
    breaks = c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = "",
    y = "% CMD+ PanEV− Anti-SA+",
    title = "bMVs–BacB Complex"
  ) +
  theme_pub_black +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.length = unit(0.5, "cm")
  )

# ===============================
# Plot: EV–Miltenyi Beads Complex
# ===============================
p2 <- ggplot(
  data,
  aes(
    Group,
    `Percentage of CMD+ PanEV+ Anti-SA+ of EV-Miltenyi beads complex`,
    fill = Group
  )
) +
  stat_summary(fun = mean, geom = "bar", width = 0.65, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1.4) +
  geom_jitter(width = 0.12, size = 3.5, shape = 21,
              color = "black", stroke = 1.3) +
  stat_pvalue_manual(
    sig_miltenyi,
    tip.length = 0.02,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  y_scale_zero +
  labs(
    x = "",
    y = "% CMD+ PanEV+ CD45+ Anti-SA+",
    title = "EV–Miltenyi Beads Complex"
  ) +
  theme_pub_black +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.length = unit(0.3, "cm")
  )

print(p1)
print(p2)

# ===============================
# Save plots
# ===============================
ggsave(
  file.path(out_dir, "CMD_PanEVneg_bMVs_significance.png"),
  p1, width = 10, height = 8, dpi = 300, bg = "white"
)

ggsave(
  file.path(out_dir, "CMD_PanEVpos_EV_significance.png"),
  p2, width = 10, height = 8, dpi = 300, bg = "white"
)

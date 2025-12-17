# ===============================
# Libraries
# ===============================
library(ggplot2)
library(dplyr)
library(ggpubr)

# ===============================
# Read data
# ===============================
data <- read.table(
  "C:/Users/ga53hil/Desktop/Bacterämia_patient_NTA.txt",
  header = TRUE,
  sep = "\t"
)

# Factor order
data$Samples <- factor(data$Samples, levels = c("HC", "BCN", "BCP"))

# ===============================
# Define group colors
# ===============================
group_colors <- c(
  HC  = "#4A4A4A",
  BCN = "#4DBBD5",
  BCP = "#E64B35"  
)

# ===============================
# Publication theme (EXTRA BLACK AXES + TEXT)
# ===============================
theme_pub_black <- theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    
    axis.title = element_text(size = 22, color = "black"),
    axis.text  = element_text(size = 20, color = "black"),
    
    axis.line  = element_line(size = 2, color = "black"),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    plot.title = element_text(
      size = 26, face = "plain", hjust = 0.5, color = "black"
    ),
    plot.subtitle = element_text(
      size = 18, face = "bold", hjust = 0.5, color = "black"
    ),
    
    legend.position = "none"
  )

# ===============================
# PLOT 1: Particle Size (nm)
# ===============================

anova_size <- aov(Particle_size_nm ~ Samples, data = data)
p_size <- summary(anova_size)[[1]][["Pr(>F)"]][1]

p1 <- ggplot(
  data,
  aes(x = Samples, y = Particle_size_nm, fill = Samples)
) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    width = 0.65,
    color = "black",
    alpha = 1
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 1.4,
    color = "black"
  ) +
  
  # Individual data points
  geom_jitter(
    width = 0.12,
    size = 3.5,
    shape = 21,
    color = "black",
    stroke = 1.3,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  labs(
    x = "Group",
    y = "Particle Size (nm)",
    title = "Particle Size",
    subtitle = paste0("One-way ANOVA, p = ", signif(p_size, 2))
  ) +
  
  theme_pub_black

print(p1)

# ===============================
# PLOT 2: CMO+ Particle Concentration
# ===============================

anova_conc <- aov(CMO_Particles_in._1_ml_serum ~ Samples, data = data)
p_conc <- summary(anova_conc)[[1]][["Pr(>F)"]][1]

p2 <- ggplot(
  data,
  aes(x = Samples, y = CMO_Particles_in._1_ml_serum, fill = Samples)
) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    width = 0.65,
    color = "black",
    alpha = 1
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 1.4,
    color = "black"
  ) +
  
  # Individual data points
  geom_jitter(
    width = 0.12,
    size = 3.5,
    shape = 21,
    color = "black",
    stroke = 1.3,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  labs(
    x = "Group",
    y = "CMO+ Particles per mL Serum",
    title = "CMO+ Particle Concentration",
    subtitle = paste0("One-way ANOVA, p = ", signif(p_conc, 2))
  ) +
  
  theme_pub_black

print(p2)

# ===============================
# Export figures (PNG, 300 dpi, white bg)
# ===============================

ggsave(
  "Particle_Size_bar_dots.png",
  p1,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "CMO_Concentration_bar_dots.png",
  p2,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

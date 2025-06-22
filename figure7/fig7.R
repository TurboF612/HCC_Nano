rm(list = ls())
gc()
# remotes::install_version("SeuratObject", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
{
  library(RColorBrewer)
  blue <- "#96C3D8"
  red <- "#F19294"
  green <- "#66C2A5"
  gray <- "#D9D9D9"
  cor_set <- c("#E7E1EF", "#C994C7", "#66C2A5", "#AEDDEE", "#80B1D3", "#3477A9", 
               "#CCEBC5", "#A4D38E", "#4A9D47", "#BADA90", "#FEE08B", "#F58135",
               "#FACB7B", "#96C3D8", "#C6DBEF", "#F5B375", "#BDA7CB",
               "#B4B1AE", "#00A43C", "#FDAE61", "#E6F598", "#FFFFB3", "#B3DE69", "#FCCDE5")
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(tidydr)
library(tidyverse)

load("8_srt_RRM2.Rda")

#### fig7a ####
fig7a <- read_tsv("fig7a.txt")

cell_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ CellType_final + RRM2_Group, 
                          data = umap_data, FUN = mean)

p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = CellType_final)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cor_set) +
  facet_wrap(~ RRM2_Group) +
  geom_text(data = cell_centers, 
            aes(x = UMAP_1, y = UMAP_2, label = CellType_final),
            size = 5, inherit.aes = FALSE) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  labs(x = "UMAP_1", y = "UMAP_2")

print(p)

ggsave("fig7a.pdf", width = 10, height = 5.5)

#### fig7b ####
fig1b <- read_tsv("fig1b.txt")
ggplot(Cell_ratio, 
       aes(Group, Freq, 
           fill = CellType, 
           stratum = CellType,
           alluvium = CellType)) +
  geom_col(width = 0.4, color = NA) +
  geom_flow(width = 0.4, alpha = 0.2, knot.pos = 0) +
  scale_fill_manual(values = cor_set) +
  theme_map() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="lightgray", size=0.5),
        axis.text.x = element_text(size = 10, vjust = 1.5, hjust = 1)) +
  RotatedAxis() +
  labs(x = "", y = "Cell Proportion") 
ggsave("fig1b.pdf", width = 6, height = 5.5)

#### fig7c-d ####
load("fig7c.Rda")
plotNhoodGraphDA(milo, milo.res.1, alpha = 0.05)
ggsave("fig7c.pdf", width = 8, height = 8)

load("fig7d.Rda")
p1 <- plotDAbeeswarm(milo.res.1, group.by = "celltype")

p1
ggsave("fig7d.pdf", width = 8, height = 8)


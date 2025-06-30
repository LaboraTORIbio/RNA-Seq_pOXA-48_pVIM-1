library(ggplot2)
library(cowplot)


# Load data
operon <- read.csv("cluster_DEGs.csv", sep="\t")

heatmap_colors <- c("gray70", "black")
operon_order <- c("ifp_pVIM-1-vs-pOXA-48", "ifp_pVIM-1-vs-PF", "ifp_pOXA-48-vs-PF",
                  "pfp_pVIM-1-vs-pOXA-48", "pfp_pVIM-1-vs-PF", "pfp_pOXA-48-vs-PF",
                  "lysR_pVIM-1-vs-pOXA-48", "lysR_pVIM-1-vs-PF", "lysR_pOXA-48-vs-PF")
operon_order <- c("ifp_pVIM-1-vs-pOXA-48", "pfp_pVIM-1-vs-pOXA-48", "lysR_pVIM-1-vs-pOXA-48",
                  "ifp_pVIM-1-vs-PF", "pfp_pVIM-1-vs-PF", "lysR_pVIM-1-vs-PF",
                  "ifp_pOXA-48-vs-PF", "pfp_pOXA-48-vs-PF", "lysR_pOXA-48-vs-PF")
strain_order <- c("K147", "K153", "K25", "KPN06", "KPN18", "K163", "K198", "K318", "KPN11", "KPN12")

# Log2FC heatmap
hmoperon <- ggplot(operon, aes(x = factor(Strain, level=strain_order), y = factor(GeneComp, level=operon_order), fill = log2FC)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#0000D5",
                       mid = "#FFFFFF",
                       high = "#D50000") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm'))

# Significance heatmap
padjoperon <- ggplot(operon, aes(x = factor(Strain, level=strain_order), y = factor(GeneComp, level=operon_order), fill = padj)) +
  geom_tile() +
  scale_fill_manual(values = heatmap_colors) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm'))



p <- plot_grid(hmoperon, padjoperon, align="v", nrow=1, ncol=2, rel_widths=c(2, 1.4), rel_heights=c(0.3, 2))
p


## Session Info

sessionInfo()

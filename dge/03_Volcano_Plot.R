suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(tidyverse)
})

deg=read.table("CTRvsKD_DESeq_ALL.txt")

data <- data.frame(gene = deg$id,
                   pvalue = -log10(deg$padj), 
                   lfc = deg$log2FoldChange,
		   abs=abs(deg$log2FoldChange))



data <- data %>%
  mutate(color = ifelse(data$lfc > 0.3 & data$pvalue > 1.3, 
                        yes = "Clock", 
                        no = ifelse(data$lfc < -0.3 & data$pvalue > 1.3, 
                                    yes = "Ctl", 
                                    no = "none")))

pdf("CLOCK_DGE_VULCANO.pdf",width=6,height=6)
# Color corresponds to fold change directionality
colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
  theme_bw(base_size = 16) + # clean up theme
  theme(legend.position = "none") + # remove legend
  xlab(expression(logFC("CLOCK" / "CTL"))) + # x-axis label
  ylab(expression(-log[10]("p-value"))) + # y-axis label
  geom_hline(yintercept = 1.3, colour = "black",linetype="dotted",size=1) + # p(0.05) = 1.3
  annotate(geom = "text", 
           label = "Down", 
           x = -1, y = 5, 
           size = 7, colour = "black",fontface = 'bold') + # add Untreated text
  annotate(geom = "text", 
           label = "Up", 
           x = 1, y = 5, 
           size = 7, colour = "black",fontface = 'bold') + # add Treated text
  scale_color_manual(values = c("Ctl" = "#E64B35", 
                                "Clock" = "#3182bd", 
                                "none" = "#636363")) # change colors
  
# Plot figure
colored

# Subset table to only show certain gene labels
sign=data[data$pvalue > 1.3,]
top_labelled <- top_n(sign, n = 10, wt = pvalue)

# Add layer of text annotation to volcano plot.
colored + geom_text_repel(data = top_labelled, 
                          mapping = aes(label = gene), 
                          size = 5,
                          fontface = 'bold', 
                          color = 'black',
                          box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.5, "lines"))
dev.off()



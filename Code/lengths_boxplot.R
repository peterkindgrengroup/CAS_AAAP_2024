library(stringr)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/")
in_beds <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/TSSseq_TES_tair10/PAS_paper"
out <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/plots/"

beds_original <-  dir(in_beds, pattern="bed", recursive = F, full.names = TRUE) 
beds_original <- beds_original[c(1,3)] #JUST PLOT CONTROLS AND PAS
#
length_subsets <- lapply (beds_original, function(x){
    table <- read.table(x, header = F, sep="\t")
    length_vector <- abs(table$V3-table$V2)
    length_vector
})
#remove path and extension
names(length_subsets) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(beds_original)) 


n.obs <- sapply(length_subsets, length)
seq.max <- seq_len(max(n.obs))
combined_table <- sapply(length_subsets, "[", i = seq.max)
combined_table <- as.data.frame(combined_table)

#statistics
ggqqplot(combined_table$PAS_noTFs)
shapiro.test(combined_table$PAS_noTFs)
#sample does not follow normal distribution, apply non-parametric test (Mann-Whitney U)
wilcox.test(combined_table$controls_all,combined_table$PAS_noTFs)$p.value 
wilcox.test(combined_table$controls_all,combined_table$PAS_TFs_all)$p.value 
wilcox.test(combined_table$PAS_noTFs,combined_table$PAS_TFs_all)$p.value 

d2.df <- reshape2::melt(as.matrix(combined_table), c("x", "y"), value.name = "z")
head(d2.df)

d2.df$y <- as.factor(d2.df$y) ## Convert the variable from a numeric to a factor variable

library(ggstatsplot)

# plot with statistical results
ggbetweenstats( # paired samples
    data = d2.df,
    x = y,
    y = z,
    type = "nonparametric", # for wilcoxon
    centrality.plotting = FALSE # remove median
)

p <- ggplot(d2.df, aes(x=y, y=z, fill=y)) + 
    geom_violin(trim=T) + labs(fill="Gene subsets")
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    ylab("Length (bp)")+
    ylim(0,6000)+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(file=paste0(out,"gene_length_control_PAS_violin.pdf"),
    useDingbats=FALSE)
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    ylab("Length (bp)")+
    ylim(0,6000)+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
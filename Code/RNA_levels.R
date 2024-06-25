library(stringr)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

#Setting working directory
setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/RNA_analysis/")
#Files with gene subsets
in_beds <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/genes_IDs"
#Output directory
out <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/RNA_analysis/plots/"

#Load RNAseq expression data
expression_control <- read.table("./RNAseq_control_peter.txt", header = F, sep="\t")
expression_control$V2 <- as.numeric(gsub("\\,", ".", expression_control$V2)) #change decimal , to . and apply numeric
expression_control <- na.omit (expression_control, cols="V2")
expression_control$V2 <- (expression_control$V2-mean(expression_control$V2))/sd(expression_control$V2) # z score

beds_original <-  dir(in_beds, recursive = F, full.names = TRUE, include.dirs = F) 
beds_original <- beds_original[ !grepl("exclude", beds_original)]

#select subsets to compare
beds_original <- beds_original[c(1,2)]


#Extract expression by gene_ID
RNA_subsets <- lapply (beds_original, function(x){
    table <- read.table(x, header = F, sep="\t")
    table <- table$V1
    subsets <- expression_control[expression_control$V1 %in% table, ]
    subsets <- subsets$V2
})
#remove path and extension
names(RNA_subsets) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(beds_original)) 

#Combine lists into dataframe
n.obs <- sapply(RNA_subsets, length)
seq.max <- seq_len(max(n.obs))
combined_table <- sapply(RNA_subsets, "[", i = seq.max)
combined_table <- as.data.frame(combined_table)

#Determine distribution
ggqqplot(combined_table$convergents)
shapiro.test(combined_table$convergents)

#sample does not follow normal distribution, apply non-parametric test (Mann-Whitney U)
wilcox.test(combined_table$controls_all,combined_table$CAS)$p.value
wilcox.test(combined_table$controls_all,combined_table$PAS_TFs_all)$p.value
wilcox.test(combined_table$PAS_TFs_all,combined_table$PAS_noTFs)$p.value

#Convert dataframe into ggplot obskect
d2.df <- reshape2::melt(as.matrix(combined_table), c("x", "y"), value.name = "z")
head(d2.df)
d2.df$y <- as.factor(d2.df$y) ## Convert the variable from a numeric to a factor variable

#plot
p <- ggplot(d2.df, aes(x=y, y=z, fill=y)) + 
    geom_violin(trim=TRUE)

#plot adjustments
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    #scale_fill_manual(values=c("#4169e1", "#663399", "#FF9900", "#CC0000")) +
    #scale_fill_manual(values=histone_colors)+
    ylim(-0.5,0.5)+
    labs(title="", x ="", y = "z-score") +
    #stat_compare_means(method = "wilcox.test", label.x = 1.3, label.y = 0.50) +
    theme_classic()

#output as pdf (repeat steps)
pdf(file=paste0(out,"RNAseq_control_PAS_zscore.pdf"),
    useDingbats=FALSE)
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    #scale_fill_manual(values=c("#4169e1", "#663399", "#FF9900", "#CC0000")) +
    #scale_fill_manual(values=histone_colors)+
    ylim(-0.5,0.5)+
    labs(title="", x ="", y = "z-score") +
    stat_compare_means(method = "wilcox.test", label.x = 1.3, label.y = 0.50) +
    theme_classic()
dev.off()
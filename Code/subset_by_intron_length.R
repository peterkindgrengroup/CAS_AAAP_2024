
#Script used to subset genes based on their 1st intron
#Allows comparison of control genes and CAS genes of similar intron length in order to find DNA motifs enriched at their intron


library(stringr)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/RNA_analysis/")
in_beds <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/fullgene_general"

expression_d7 <- read.table("./RNA_day7_cellculture.txt", header = F, sep="\t")
expression_d7$V2 <- as.numeric(gsub("\\,", ".", expression_d7$V2)) #change decimal , to . and apply numeric
expression_d7 <- na.omit (expression_d7, cols="V2")
expression_d7$V2 <- (expression_d7$V2-mean(expression_d7$V2))/sd(expression_d7$V2)

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/intron_analysis")
first_introns<-read.csv(file = "first_introns_araport.txt",sep = "\t", header = F)



#automate intron length count
gene_IDs <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/genes_IDs"

beds_original <-  dir(gene_IDs, recursive = F, full.names = TRUE) 
beds_original <- beds_original[ !grepl("exclude", beds_original) ]

#keep GENE_ID and intron length value of each category, add expression data
subsets_1st_intron_RNA <- lapply (beds_original, function(x){
    subsets <- scan(file = x, character(), quote = "")
    subsets_introns <- first_introns[first_introns$V1 %in% subsets, ]
    subsets_RNA <- expression_d7[expression_d7$V1 %in% subsets, ]
    joined_subsets <- merge(subsets_introns, subsets_RNA, by.x = "V1", by.y = "V1", all.x = TRUE, all.y = FALSE)
    joined_subsets
})

#remove path and extension
names(subsets_1st_intron_RNA) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(beds_original)) 

for(x in 1:length(subsets_1st_intron_RNA)) {
    elemento <- subsets_1st_intron_RNA[[x]]
    nombre <- names(subsets_1st_intron_RNA)[x]
    
    # Write table of gene IDs
    write.table(elemento, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/intron_analysis/split_by_1stintron_RNA/",nombre,"_1st_intron_RNA.txt"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
    
}

#After filter for 1st length intron (300,3000) and RNA comparable levels
subset_IDs <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/intron_analysis/subset_by_intron_RNA"
subset_files <-  dir(subset_IDs, recursive = F, full.names = TRUE) 
subset_files <- subset_files[ !grepl("exclude", subset_files) ]

controls_TSS_TAIR <- read.table(file = "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/TSSseq_TES_tair10/controls_all.bed", sep = "\t", quote = "")
controls_subseted <- read.table(file = subset_files[2], sep = "\t", quote = "")
controls_TSS_TAIR <- controls_TSS_TAIR[controls_TSS_TAIR$V5 %in% controls_subseted$V1, ]
write.table(controls_TSS_TAIR, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/subsets_intronsControl/controls_intronControl.bed"), sep ="\t", row.names = F, col.names = F, quote=FALSE )


CAS_TSS_TAIR <- read.table(file = "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/TSSseq_TES_tair10/convergents.bed", sep = "\t", quote = "")
CAS_subseted <- read.table(file = subset_files[4], sep = "\t", quote = "")
CAS_TSS_TAIR <- CAS_TSS_TAIR[CAS_TSS_TAIR$V5 %in% CAS_subseted$V1, ]
write.table(CAS_TSS_TAIR, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/beds/subsets_intronsControl/convergents_intronControl.bed"), sep ="\t", row.names = F, col.names = F, quote=FALSE )




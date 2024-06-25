library(stringr)
library(data.table)
library(dplyr)
length_subsets <- data.frame(matrix(ncol = 2, nrow = 0))

# setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/")
# out <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/"

setwd("~/MEGA/Work/Shiv et al 2022/RNA_analysis/input_subsets/genes_IDs/")
out="~/MEGA/Work/Shiv et al 2022/RNA_analysis/input_subsets/genes_IDs/"

TSSseq <- read.table(file = "~/MEGA/Work/Shiv et al 2022/RNA_analysis/input_subsets/TSSseq_arabidopsis.txt" , header = F, sep="\t")
TSSseq <- TSSseq[,c(1:6)]

TAIR10 <- read.table(file = "~/MEGA/Work/Shiv et al 2022/RNA_analysis/input_subsets/TAIR10_genes_transposons.bed12" , header = F, sep="\t")
TAIR10 <- TAIR10[,c(1:6)]
TAIR10$V4 <- str_sub(TAIR10$V4, end=-3)

name_of_subset <- "PAS_TFs_all"
subsets <- scan(file = paste0(name_of_subset), character(), quote = "")
subsets_TAIR <- subsets #save original for later
# subset from TSSseq
subsets <- TSSseq[TSSseq$V5 %in% subsets, ]
#remove duplicates created by lnRNA features,etc
subsets <- subsets[!duplicated(subsets$V5), ]
colnames(subsets)<-c("V1","V2","V3","V4","geneID","V6")
#NEED TO ADD Chr for Homer:
subsets$V1 <- sub("^", "Chr", subsets$V1 )


# subset from TAIR10
subsets_TAIR <- TAIR10[TAIR10$V4 %in% subsets_TAIR, ]
#remove duplicates created by lnRNA features,etc
subsets_TAIR <- subsets_TAIR[!duplicated(subsets_TAIR$V4), ]
colnames(subsets_TAIR)<-c("V1","V2","V3","geneID","V5","V6")

#TSSseq is limiting as not all genes are included
TSSseq_ttsTAIR <- merge(x = subsets, y = subsets_TAIR, by = "geneID", all.x = TRUE)

### considered USING CENTER OF TSSseq as START, no, it fails on - strand ones
## CONSIDER using left for positives and right for negatives
colnames(TSSseq_ttsTAIR) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11")

positives <- filter(TSSseq_ttsTAIR, V6 == "+")
positives <- positives[,c(2,3,9,5,1,6)]

negatives <- filter(TSSseq_ttsTAIR, V6 == "-")
negatives <- negatives[,c(2,8,4,5,1,6)]
colnames(negatives)<-colnames(positives)

TSSseq_ttsTAIR <- rbindlist(list(positives, negatives))
TSSseq_ttsTAIR <- TSSseq_ttsTAIR[order(V1,V2),]
write.table(TSSseq_ttsTAIR, file = paste0(out,"2024_coordinates/004_",name_of_subset,".bed"), sep ="\t", row.names = F, col.names = F, quote=FALSE )


#########


mean_length <- mean(abs(TSSseq_ttsTAIR$V3-TSSseq_ttsTAIR$V9))
std <- function(x) sd(x)/sqrt(length(x))
se_length <- std(abs(TSSseq_ttsTAIR$V3-TSSseq_ttsTAIR$V9))
columns_toappend <- as.matrix(t(c(mean_length,se_length,name_of_subset)))
# 
# in case to calculate the number of observations
# mean_length <- length(TSSseq_ttsTAIR$V3)
# columns_toappend <- as.matrix(t(c(mean_length,name_of_subset)))

combined_lengths <- rbind(combined_lengths,columns_toappend)

#rename files
#rename -n 's/trimmomatic.sorted_nolowup_extended_50bp.bw//' *.gz


#AT some point plot lenghts of each category!
#subsets$length <- as.integer(abs(subsets$V3 - subsets$V2))


#USE BED12 in /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/TAIR10_GFF3_genes_transposons.gff.bed12
#BED6 format Chr start end name score strand
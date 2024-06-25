setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/intron_analysis")

library(dplyr)
library(stringr)
library(GenomicFeatures) 
library(rtracklayer)

#raw data from Araport11 with introns
txdb <- makeTxDbFromGFF(file = "Araport11_GFF3_genes_transposons_introns.201606.gff", format = "gff") 
introns <- intronsByTranscript(txdb, use.names=TRUE)
geneid = mapIds(txdb, names(introns), "GENEID", "TXNAME")

intron_GRange <- unlist(introns)

negative_strand <- intron_GRange[strand(intron_GRange) == "-"]
positive_strand <- intron_GRange[strand(intron_GRange) == "+"]

positives <- data.frame(geneID=names(positive_strand),seqnames=seqnames(positive_strand),
                        starts=start(positive_strand)-1,
                        ends=end(positive_strand),
                        strands=strand(positive_strand))
#same as before
positives <- distinct(positives, geneID,.keep_all= TRUE)
positives <- distinct(positives, seqnames, starts, ends, .keep_all= TRUE)


negatives <- data.frame(geneID=names(negative_strand),seqnames=seqnames(negative_strand),
                 starts=start(negative_strand)-1,
                 ends=end(negative_strand),
                 strands=strand(negative_strand),
                 sort=seq(length.out=length(negative_strand)))
#reverse list to take last exon
negatives <- negatives[order(negatives[,6], decreasing = TRUE),]
#same as before
negatives <- distinct(negatives, geneID,.keep_all= TRUE)
negatives <- distinct(negatives, seqnames, starts, ends, .keep_all= TRUE)

#reorder
negatives <- negatives[order(negatives[,6], decreasing = FALSE),]
#remove last column so it can be rbinded
negatives <- negatives[,c(1:5)]

#join 2 populations
first_introns <- rbind(positives,negatives)

#calculate 
first_introns$length <- abs(first_introns$ends-first_introns$starts)
first_introns$geneID <- str_sub(first_introns$geneID, end=-3)
write.table(first_introns,file = "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/intron_analysis/first_introns_araport.txt",sep = "\t",row.names = F, col.names = F, quote = F )


#automate intron length count
gene_IDs <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/input_subsets/genes_IDs"
out <- "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/2022_Peter/plots/"

beds_original <-  dir(gene_IDs, recursive = F, full.names = TRUE) 
beds_original <- beds_original[ !grepl("exclude", beds_original) ]

#keep just the exon length value of each category
length_1st_intron <- lapply (beds_original, function(x){
    subsets <- scan(file = x, character(), quote = "")
    subsets <- first_introns[first_introns$geneID %in% subsets, ]
    length_vector <- subsets$length
    length_vector
    
})
#remove path and extension
names(length_1st_intron) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(beds_original)) 

n.obs <- sapply(length_1st_intron, length)
seq.max <- seq_len(max(n.obs))
combined_table <- sapply(length_1st_intron, "[", i = seq.max)
combined_table <- as.data.frame(combined_table)

d2.df <- reshape2::melt(as.matrix(combined_table), c("x", "y"), value.name = "z")
head(d2.df)

d2.df$y <- as.factor(d2.df$y) ## Convert the variable from a numeric to a factor variable

p <- ggplot(d2.df, aes(x=y, y=z, fill=y)) + 
    geom_violin(trim=T) + labs(fill="Gene subsets")
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    ylab("Length (bp)")+
    ylim(0,750)+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(file=paste0(out,"first_intron_length.pdf"),
    useDingbats=FALSE)
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    ylab("Length (bp)")+
    ylim(0,750)+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()





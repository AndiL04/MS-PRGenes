###iRIGS
####Step 0
###supporting file processing



###--Step 0.1-Generating differential expression information for MS---------###
### Loading packages###
library(DESeq2)

###Working address
setwd("C:/Users/andya/OneDrive/Desktop/MS Bayes&2SMR/iRIGS")

####Read Data###
ms_gene_quant <- read.delim("./supporting_files/GSE111972_norm_data.txt",
                            row.names = 1, stringsAsFactors = T)
ms_gene_quant <- ms_gene_quant[,grepl("WM",colnames(ms_gene_quant))]

###Round count###
ms_rounded_count <- floor(ms_gene_quant + 0.5)

##Gnerating phenotype information###
phenotype <- data.frame(row.names = colnames(ms_rounded_count),
                        condition = colnames(ms_rounded_count))

phenotype[grepl("CON",row.names(phenotype)),] <- "Control"
phenotype[!grepl("CON",row.names(phenotype)),] <- "MS"


###Differential Analysis with DESeq2###
dds <- DESeqDataSetFromMatrix(countData = ms_rounded_count,
                              colData = phenotype,
                              design= ~ condition)
##Set ref level with DESeq2###
dds$Condition <- relevel(dds$condition, ref = 'Control')

###Proceeding###
dds <- DESeq(dds)
resultsNames(dds)

###Output result###
MS_differential <- results(dds, name = "condition_MS_vs_Control")
MS_differential_ordered <- MS_differential[order(MS_differential$padj),]
write.csv(MS_differential_ordered, file = "./supporting_files/MS_DE.csv") 
###"MS_DE.CSV" would be processed in iRIGS.R




###--Step 0.2-Generating differential methylation information for iRGIS-----###

### Loading Data###
methy.z <- read.delim("./supporting_files/MethyNodeWeights_MSNAWM.txt", header = F)

methy.z$p_value <-  2*pnorm(methy.z$V2, lower.tail = F)


colnames(methy.z) <- c("gene","z_score","p_value")
write.csv(methy.z, "./supporting_files/methylation.ms.csv")




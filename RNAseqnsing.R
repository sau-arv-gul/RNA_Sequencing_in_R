if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RUVSeq")

require("RUVSeq")
s
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

require("DESeq2")
s
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("zebrafishRNASeq")

require("zebrafishRNASeq")
s

library(RUVSeq)
library(zebrafishRNASeq)
data(zfGenes)
head(zfGenes)
tail(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)
set1 <- RUVg(set, spikes, k=1)
pData(set1)

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)

design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
res

output <- as.data.frame(res)

final <- data.frame(row.names(output))
final[,c(2,3)] <- output[,c(2,5)]

names(final)[1] <- "Gene"

with(final, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Saurabh Kumar_2020541_DEGs"))

with(subset(final, pvalue < .05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))
 
#lets store the whole data into DATA.csv file
df <- data.frame(res)
df
write.csv(df,"DATA.csv") 
library(dplyr)
sample1 = filter(df,pvalue<0.05)  # sample1 contain all genes with p<0.05
write.csv(sample1,"data_less.csv")  #storing sample1 in csv file
s2 = filter(df,log2FoldChange > 0) #sample for upregulated2
s2 <- s2 %>% arrange(desc(log2FoldChange))
s2 <- head(s2, 50)
write.csv(s2,"TOP_UP_REGULATED.csv")            # top upregulated 
s3 <- filter(df,log2FoldChange < 0)
s3 <-  s3 %>% arrange(log2FoldChange)
s3 <- head(s3, 50)
write.csv(s3,"TOP_DOWN_REGULATED.csv")



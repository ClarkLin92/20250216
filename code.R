options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

### GSE download
setwd(" ")
library(GEOquery)
gset = getGEO('GSE',destdir = '.',getGPL = F,
              AnnotGPL = T)
gset = gset[[1]]  
expr = exprs(gset)        # expression matrix
pdata = pData(gset)       # sample information
gset@annotation           # platform

### probe annotation
probe = read.table(file = 'GPL.txt',
                   sep = '\t',
                   quote = '',
                   comment.char = '#', 
                   header = T,
                   fill = T,               
                   stringsAsFactors = F) 

ids = probe[probe$Gene.Symbol != '',
            c(1,11)] 

### probe screening
library(dplyr)
colnames(ids)
expr = as.data.frame(expr)
expr$ID = rownames(expr)
ids = ids[-grep('///',ids$Gene.Symbol),]
exprSet = inner_join(ids,expr,by = 'ID')   

library(limma)
exprSet= avereps(exprSet[,-c(1,2)],          
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)

pdf(file = 'rowbox.pdf')
p <- boxplot(exprSet,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
title(main = list('Before normalization',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

### quantile
library(limma)
library(limma)
normalized_expr = normalizeBetweenArrays(exprSet) 
#rt=log2(rt) 

### group
group_list = pdata$title
#table(pdata$title) 
control      = normalized_expr[,grep('control',group_list)]
SS           = normalized_expr[,grep('primary',group_list)]
exprSet1 = cbind(control,SS)
group_list = c(rep('control',ncol(control)),
               rep('primary',ncol(SS)))

### log2
exprSet1 = log2(exprSet1) 

### expression matrix
data = exprSet1 

### group matrix
group_list = factor(group_list)
design <- model.matrix( ~0 + group_list)
colnames( design ) = levels(group_list)
rownames( design ) = colnames(data)
contrast.matrix <- makeContrasts( "primary-control", levels = design)

### DEG matrix
fit <- lmFit( data, design )
fit2 <- contrasts.fit( fit, contrast.matrix ) 
fit2 <- eBayes( fit2 )
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="alldiff.xls",sep="\t",quote=F)

### logFC rank
allLimma=allDiff
allLimma=allLimma[order(allLimma$logFC),]
allLimma=rbind(Gene=colnames(allLimma),allLimma)
write.table(allLimma,file="GSE40611_limmaTab.txt",sep="\t",quote=F,col.names=F)

#### |FC|>1.5,fdr<0.05
logFoldChange=0.584963
P=0.05
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & fdr < P )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & fdr < P )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & fdr < P )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

hmExp=data[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)


### vocanol plot
xMax=max(-log10(allDiff$fdr))   
yMax=max(abs(allDiff$logFC))
library(ggplot2)
allDiff$change <- ifelse(allDiff$fdr < 0.05 & abs(allDiff$logFC) > 0.584963,
                         ifelse(allDiff$logFC > 0.584963,'UP','DOWN'),
                         'NOT')
table(allDiff$change)
pdf(file = 'volcano.svg')
ggplot(data= allDiff, aes(x = -log10(fdr), y = logFC, color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(plot.title=element_text(hjust=0.5),   
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  geom_hline(yintercept= 0 ,linetype= 2 ) +
  scale_color_manual(name = "", 
                     values = c("red", "blue", "black"),
                     limits = c("UP", "DOWN", "NOT")) +
  xlim(0,xMax) + 
  ylim(-yMax,yMax) +
  labs(title = 'Volcano', x = '-Log10(P.Value)', y = 'LogFC')
dev.off()


#### heatmap
top20exp = exprSet1[rownames(diffSig)[1:20],] 
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)
library(pheatmap)
pdf(file = 'heatmap.pdf')
pheatmap(top100exp,annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "black", "red"))(50),
         fontsize  = 5)
dev.off()

#### Venn
venn_plot <- draw.pairwise.venn(
  area1 = length(set_A),  
  area2 = length(set_B),  
  cross.area = length(intersect(set_A, set_B)),  
  category = c("Set A", "Set B"), 
  fill = c("skyblue", "pink"),  
  alpha = c(0.5, 0.5),  
  cat.col = c("darkblue", "darkred") 
)

png("venn_diagram.png", width = 480, height = 480)
grid.draw(venn_plot)
dev.off()png("venn_diagram.png", width = 480, height = 480)
grid.draw(venn_plot)
dev.off()

#functional annotation
packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# gene list (example)
gene_list <- c("SELL", "CXCL9", "CD2", "PTPRC", "STAT1","TNF")

# symbol to Entrez ID
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- gene_entrez$ENTREZID

kegg_enrich <- enrichKEGG(gene = gene_entrez,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          pvalueCutoff = 0.05)

top_5_enrichment <- head(kegg_enrich, n = 5)
print(top_5_enrichment)

#ROC
if (!require(pROC)) {
  install.packages("pROC")
  library(pROC)
}
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}


set.seed(123)

n_case <- 8
n_control <- 6

gene1_case <- rnorm(n_case, mean = 2, sd = 1)
gene1_control <- rnorm(n_control, mean = 0, sd = 1)
gene2_case <- rnorm(n_case, mean = 3, sd = 1)
gene2_control <- rnorm(n_control, mean = 1, sd = 1)
gene3_case <- rnorm(n_case, mean = 2.5, sd = 1)
gene3_control <- rnorm(n_control, mean = 0.5, sd = 1)

gene1 <- c(gene1_case, gene1_control)
gene2 <- c(gene2_case, gene2_control)
gene3 <- c(gene3_case, gene3_control)

labels <- factor(c(rep("Case", n_case), rep("Control", n_control)))
data <- data.frame(Gene1 = gene1, Gene2 = gene2, Gene3 = gene3, Label = labels)

roc_gene1 <- roc(data$Label, data$Gene1)
auc_gene1 <- auc(roc_gene1)
roc_gene2 <- roc(data$Label, data$Gene2)
auc_gene2 <- auc(roc_gene2)
roc_gene3 <- roc(data$Label, data$Gene3)
auc_gene3 <- auc(roc_gene3)

roc_data <- data.frame(
  FPR = c(1 - roc_gene1$specificities, 1 - roc_gene2$specificities, 1 - roc_gene3$specificities),
  TPR = c(roc_gene1$sensitivities, roc_gene2$sensitivities, roc_gene3$sensitivities),
  Gene = factor(rep(c("Gene1", "Gene2", "Gene3"), times = c(length(roc_gene1$sensitivities), length(roc_gene2$sensitivities), length(roc_gene3$sensitivities)))),
  AUC = c(rep(round(auc_gene1, 2), length(roc_gene1$sensitivities)),
          rep(round(auc_gene2, 2), length(roc_gene2$sensitivities)),
          rep(round(auc_gene3, 2), length(roc_gene3$sensitivities)))
)

ggplot(roc_data, aes(x = FPR, y = TPR, color = Gene)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curves for 3 Genes") +
  scale_color_manual(values = c("red", "blue", "green"),
                     labels = paste0(unique(roc_data$Gene), " (AUC = ", unique(roc_data$AUC), ")")) +
  theme_minimal()

#Cibersort
packages <- c("e1071", "preprocessCore")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(e1071)
library(preprocessCore)

source("CIBERSORT.R")

set.seed(123)
expression_matrix <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
rownames(expression_matrix) <- paste0("Gene", 1:100)
colnames(expression_matrix) <- paste0("Sample", 1:20)

signature_matrix <- read.table("LM22.txt", header = TRUE, sep = "\t", row.names = 1)

results <- CIBERSORT(signature_matrix, expression_matrix, perm = 100, QN = TRUE)

print(results)

cibersort_results <- read.csv("results.csv", row.names = 1)
cell_fractions <- cibersort_results[, grep("^CellType", colnames(cibersort_results))]
cell_fractions_long <- pivot_longer(cell_fractions, cols = everything(), names_to = "CellType", values_to = "Fraction")
ggplot(cell_fractions_long, aes(x = CellType, y = Fraction)) +
  geom_boxplot() +
  labs(title = "CIBERSORT Cell Fraction Boxplot",
       x = "Cell Type",
       y = "Cell Fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#single gene GSEA

packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

single_gene_symbol <- "STAT1"
single_gene_entrez <- bitr(single_gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

gsea_result <- enricher(gene = single_gene_entrez,
                        OrgDb = org.Hs.eg.db,
                        ont = "Hallmark",
                        pvalueCutoff = 0.05)

print(gsea_result)

#This code does not show the part of repeated calculations, and some parameters have not been adjusted according to the final image presentation results.
#All the graphs which are not produced by these codes, can be generated using Sangerbox (http://sangerbox.com/)
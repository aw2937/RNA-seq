library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.At.tair.db)
library(RColorBrewer)
library(ape)
library(biomaRt)
library(dplyr)
library(ReactomePA)
library(clusterProfiler)

# Set working directory to source file location
setwd("/Volumes/Transcend/p53/RNA_seq")
#############
# Read data
raw_counts <- read.delim("counts.txt", stringsAsFactors = FALSE,comment.char = "#")
rownames(raw_counts) <- raw_counts[,1]
counts <- raw_counts[,-(1:6)]

# Read sample information
sampleinfo <- read.delim("sample_info.txt")
colnames(counts)<-c(sampleinfo$Sample)

sampleinfo$Treatment<-factor(sampleinfo$Treatment)
sampleinfo<-sampleinfo[order(sampleinfo$Sample,decreasing = T),]
counts<-counts[,order(colnames(counts),decreasing = T)]

#############
#############
# Count processing via CPM for filtering
CPM <- cpm(counts)
# Filter highly expressed genes
thresh <- CPM > 0.5
high <- rowSums(thresh) >= 2
counts_high <- counts[high,]
# Store filtered count data in DGEList object
y <- DGEList(counts_high)
# Plot library sizes
pdf("images/library_sizes.pdf",6,4) 
par(mar = c(8, 5, 5, 5))
barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Library sizes")
dev.off() 
# Get log2 counts per million
log_counts <- cpm(y,log=TRUE)
#############
# MDSplot of principle components analysis
col_cell <- c("red","green","blue")[sampleinfo$Treatment]
pdf("images/MDS.pdf",6,4) 
plotMDS(y,col=col_cell,pch=16)
legend("bottomright",fill=c("red","green","blue"),legend=levels(sampleinfo$Treatment))
title("MDS of Treatment")
dev.off() 
#############
# Heatmap for the 100 most variable genes
var_genes <- apply(log_counts, 1, var)
var_select <- names(sort(var_genes, decreasing=TRUE))[1:100]
var_high <- log_counts[var_select,]
pal <- brewer.pal(11,"RdYlBu")
more_col <- colorRampPalette(pal)
col_cell <- c("red","green","blue")[sampleinfo$Treatment]
pdf("images/heatmap_500.pdf",10,10) 
par(mar = c(5, 2, 2, 2))
heatmap.2(var_high,col=rev(more_col(50)),trace="none", main="Top 100 most variable genes",ColSideColors=col_cell,scale="row",srtCol=30)
dev.off() 
#############
# Normalisation for composition bias
y <- calcNormFactors(y)
# Plot changes
pdf("images/Normalization.pdf",12,12) 
par(mfrow=c(2,2))
plotMD(log_counts,column = 2)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")
plotMD(log_counts,column = 8)
abline(h=0,col="grey")
plotMD(y,column = 8)
abline(h=0,col="grey")
dev.off() 
#############
# Read and create groups
group <- sampleinfo$Treatment
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#############
#Use voom for normalization
voom <- voom(y,design,plot = FALSE)
# Plot distributions of samples before and after normalization
pdf("images/Normalization_comparison.pdf",12,6) 
par(mfrow=c(1,2))
par(mar = c(8, 5, 5, 5))
boxplot(log_counts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(log_counts),col="blue")
boxplot(voom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(voom$E),col="blue")
dev.off() 
#############
# Differential Expression: fit linear model for each gene
fit <- lmFit(voom)
# Create contrast matrix and apply to fit object
cont_matrix <- makeContrasts(ctrlVsAZA= AZA-cntrl , ctrlVsDAC= DAC-cntrl,levels=design)
fit_contrast <- contrasts.fit(fit, cont_matrix)
fit_contrast <- eBayes(fit_contrast)
sum_fit <- decideTests(fit_contrast)
summary(sum_fit)
#############
# Plot MD and volcano with highlight of 20 most DE genes
pdf("images/MD_volcano.pdf",12,12) 
par(mfrow=c(2,2))
plotMD(fit_contrast,coef=1,status=sum_fit[,"ctrlVsAZA"], values = c(-1, 1))
volcanoplot(fit_contrast,coef=1,highlight=20,names=rownames(fit_contrast))
title("20 most DE genes (ctrlVsAZA)")
plotMD(fit_contrast,coef=2,status=sum_fit[,"ctrlVsDAC"], values = c(-1, 1))
volcanoplot(fit_contrast,coef=2,highlight=20,names=rownames(fit_contrast))
title("20 most DE genes (ctrlVsDAC)")
dev.off() 
##
# Add ENTREZID to annotation

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(fit_contrast)
df<-data.frame(entrezgene_id=genes)
G_list <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=genes,mart= mart)
G_list$entrezgene_id<-as.character(G_list$entrezgene_id)
idmap <- left_join(x = df, y = G_list, by = "entrezgene_id")
idmap<-idmap[!duplicated(idmap$entrezgene_id), ]
fit_contrast$genes <- idmap$hgnc_symbol

#############
# Plot normalised log2 expression of individual genes (eg. microtubule-binding proteins CC1/MAP65-1)
pdf("Figures/R/gene_expression.pdf",6,6) 
par(mfrow=c(2,2))
nice.col <- brewer.pal(6,name="Dark2")
stripchart(voom$E["AT5G52310",]~factor(group, levels = c("ctrl","cold_1h", "cold_24h")),vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="LTI78")
stripchart(voom$E["AT5G25130",]~factor(group, levels = c("ctrl","cold_1h", "cold_24h")),vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Cytochrome P450")
stripchart(voom$E["AT1G45688",]~factor(group, levels = c("ctrl","cold_1h", "cold_24h")),vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="CC1")
stripchart(voom$E["AT1G17620",]~factor(group, levels = c("ctrl","cold_1h", "cold_24h")),vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="AT1G17620")
dev.off() 
#############
# Test against treshold
fit_treat <- treat(fit_contrast,lfc=1) #logFC of 1 = fold change of 2 or 0.5
res_treat <- decideTests(fit_treat)
summary(res_treat)

all_fit<-topTable(fit_treat,coef=1,sort.by="p",n=Inf)
all_fit$col<-rep("no_sig")
sig<-topTable(fit_treat,coef=1,sort.by="p",p.value = 0.05,lfc=1,n=Inf)
all_fit[rownames(all_fit) %in% rownames(sig),]$col<-"sig"
aza_plot <- ggplot(data=all_fit, aes(x=logFC, y=-log10(adj.P.Val),col=col))  + geom_point()  + xlim(-5,10) + ylim(0, 3) + theme_minimal() + scale_color_manual(values = c("no_sig" = "grey", "sig" = "#F7C224"))

all_fit<-topTable(fit_treat,coef=2,sort.by="p",n=Inf)
all_fit$col<-rep("no_sig")
sig<-topTable(fit_treat,coef=2,sort.by="p",p.value = 0.05,lfc=1,n=Inf)
all_fit[rownames(all_fit) %in% rownames(sig),]$col<-"sig"
dac_plot <- ggplot(data=all_fit, aes(x=logFC, y=-log10(adj.P.Val),col=col)) + geom_point()+ xlim(-5,10) + ylim(0, 3)   + theme_minimal() + scale_color_manual(values = c("no_sig" = "grey", "sig" = "#109CC8"))

require(gridExtra)
grid.arrange(aza_plot, dac_plot, ncol=2)


# Replot after threshold
pdf("images/MD_after_threshold.pdf",8,8) 
par(mfrow=c(2,2))
plotMD(fit_contrast,coef=1,status=sum_fit[,"ctrlVsAZA"], values = c(-1,1))
abline(h=0,col="grey")
plotMD(fit_contrast,coef=2,status=sum_fit[,"ctrlVsDAC"], values = c(-1,1))
abline(h=0,col="grey")
plotMD(fit_treat,coef=1,status=res_treat[,"ctrlVsAZA"], values=c(-1,1))
abline(h=0,col="grey")
plotMD(fit_treat,coef=2,status=res_treat[,"ctrlVsDAC"], values=c(-1,1))
abline(h=0,col="grey")
dev.off()

pdf("images/MD_volcano.pdf",12,12) 
par(mfrow=c(2,1))
volcanoplot(fit_treat,coef=1,highlight=111,names=rownames(fit_treat),xlim=c(-5,10),ylim=c(0,10))
title("20 most DE genes (ctrlVsAZA)")
volcanoplot(fit_treat,coef=2,highlight=13,names=rownames(fit_treat),xlim=c(-5,10),ylim=c(0,10))
title("20 most DE genes (ctrlVsDAC)")
dev.off() 

#############
# Interactive MD & Volcano HTML plots
glXYPlot(x=fit_contrast$coefficients[,1], y=fit_contrast$lods[,1],
         xlab="logFC", ylab="B", main="ctrlVsAZA",
         counts=y$counts, group=factor(group, levels = c("cntrl","AZA", "DAC")), status=sum_fit[,1],
         side.main="geneID", folder="volcano")
glMDPlot(fit_treat, coef=1, counts=y$counts, group=factor(group, levels = c("cntrl","AZA", "DAC")),
         status=res_treat, side.main="geneID", main="ctrlVsAZA",
         folder="AZA")
glMDPlot(fit_treat, coef=2, counts=y$counts, group=factor(group, levels = c("cntrl","AZA", "DAC")),
         status=res_treat, side.main="geneID", main="ctrlVsDAC",
         folder="DAC")

#Export
DEgenes<-list()
#DE<-topTable(fit_contrast,coef=1,sort.by="p",p.value = 0.05,lfc=1,n=Inf)
DE<-topTable(fit_treat,coef=1,sort.by="p",p.value = 0.05,lfc=1,n=Inf)
DEgenes[["AZA_up"]]<-rownames(DE[DE$logFC > 0,])
DEgenes[["AZA_down"]]<-rownames(DE[DE$logFC < 0,])

DE<-topTable(fit_treat,coef=2,sort.by="p",p.value = 0.05,lfc=1,n=Inf)
DEgenes[["DAC_up"]]<-rownames(DE[DE$logFC > 0,])
DEgenes[["DAC_down"]]<-rownames(DE[DE$logFC < 0,])



pathway <- lapply(DEgenes,enrichPathway,organism="human")
pathway<-pathway[lapply(pathway,function(i) length(i$ID))>0]

KEGG<-lapply(DEgenes,enrichKEGG,organism="hsa")
KEGG<-KEGG[lapply(KEGG,function(i) length(i$ID))>0]

GO<-lapply(DEgenes,enrichGO,OrgDb="org.Hs.eg.db")
GO<-GO[lapply(GO,function(i) length(i$ID))>0]

lapply(seq_along(pathway),function(i) dotplot(pathway[[i]],showCategory = 30,title=paste("Pathway enrichment of ",(names(pathway)[[i]]))))
lapply(seq_along(KEGG),function(i) dotplot(KEGG[[i]],showCategory = 10,title=paste("KEGG enrichment of ",(names(KEGG)[[i]]))))
lapply(seq_along(GO),function(i) dotplot(GO[[i]],showCategory = 10,title=paste("GO enrichment of ",(names(GO)[[i]]))))

venn(DEgenes)

all_up<-c(AZA_DAC_up,AZA_only_up)

AZA_DAC_up<-intersect(DEgenes$AZA_up,DEgenes$DAC_up)
AZA_only_up<-setdiff(DEgenes$AZA_up,DEgenes$DAC_up)
DAC_only_up<-setdiff(DEgenes$DAC_up,DEgenes$AZA_up)
all_up<-c(AZA_DAC_up,AZA_only_up)


commonpath<-enrichPathway(AZA_DAC_up, organism="human")
dotplot(commonpath,showCategory = 10)
commonKEGG<-enrichKEGG(AZA_DAC_up, organism="human")
dotplot(commonKEGG,showCategory = 10)
commonGO<-enrichGO(AZA_DAC_up, OrgDb="org.Hs.eg.db")
dotplot(commonGO,showCategory = 10)


commonpath<-enrichPathway(all_up, organism="human")
dotplot(commonpath,showCategory = 10)

commonpath<-enrichPathway(DAC_only_up, organism="human")
dotplot(commonpath,showCategory = 10)


names <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","hgnc_symbol"),values=AZA_DAC_up,mart= mart)


grid.newpage()
draw.pairwise.venn(area1 = 111, area2 = 13, cross.area = 13, category = c("AZA", 
                                                                       "DAC"), lty = rep("blank", 2), fill = c("#F7C224", "#109CC8"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 180), euler.d = TRUE, sep.dist = 0.03, 
                   rotation.degree = 45)

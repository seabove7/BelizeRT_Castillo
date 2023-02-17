#####Ssid reciprocal transplant RNAseq analysis
###host data
#setwd("~/Dropbox/UNC/RT_Paper/2022_Final_RT/GE")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("vegan")
library("tidyverse")
library("pheatmap")
library("VennDiagram")
library("adegenet") 
library("WGCNA")

# sessionInfo()
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] WGCNA_1.61                 fastcluster_1.1.24         dynamicTreeCut_1.63-1     
# [4] adegenet_2.1.0             ade4_1.7-8                 VennDiagram_1.6.18        
# [7] futile.logger_1.4.3        pheatmap_1.0.8             forcats_0.5.1             
# [10] stringr_1.4.0              purrr_0.3.4                readr_1.4.0               
# [13] tidyr_1.1.3                tibble_3.1.2               tidyverse_1.3.1           
# [16] vegan_2.5-4                lattice_0.20-35            permute_0.9-4             
# [19] dplyr_1.0.7                ggplot2_3.3.5              DESeq2_1.16.1             
# [22] SummarizedExperiment_1.6.5 DelayedArray_0.2.7         matrixStats_0.52.2        
# [25] Biobase_2.36.2             GenomicRanges_1.28.6       GenomeInfoDb_1.12.3       
# [28] IRanges_2.10.5             S4Vectors_0.14.7           BiocGenerics_0.22.1       
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1            backports_1.1.2         Hmisc_4.1-0             plyr_1.8.4             
# [5] igraph_1.1.2            sp_1.2-7                splines_3.4.2           BiocParallel_1.10.1    
# [9] robust_0.4-18           digest_0.6.18           foreach_1.4.3           htmltools_0.3.6        
# [13] GO.db_3.4.1             gdata_2.18.0            fansi_0.4.0             magrittr_2.0.1         
# [17] checkmate_1.8.5         memoise_1.1.0           fit.models_0.5-14       doParallel_1.0.11      
# [21] cluster_2.0.6           annotate_1.54.0         modelr_0.1.8            gmodels_2.16.2         
# [25] colorspace_1.4-0        rrcov_1.4-3             blob_1.2.1              rvest_1.0.0            
# [29] haven_2.4.1             xfun_0.24               crayon_1.4.1            RCurl_1.95-4.8         
# [33] jsonlite_1.7.2          genefilter_1.58.1       impute_1.50.1           survival_2.41-3        
# [37] iterators_1.0.8         ape_5.5                 glue_1.4.2              gtable_0.2.0           
# [41] zlibbioc_1.22.0         XVector_0.16.0          seqinr_3.4-5            DEoptimR_1.0-8         
# [45] scales_1.0.0            mvtnorm_1.0-7           futile.options_1.0.0    DBI_1.1.1              
# [49] Rcpp_1.0.0              xtable_1.8-2            spData_0.2.6.9          htmlTable_1.11.1       
# [53] foreign_0.8-69          bit_1.1-12              spdep_0.7-4             preprocessCore_1.38.1  
# [57] Formula_1.2-2           htmlwidgets_1.3         httr_1.4.2              RColorBrewer_1.1-2     
# [61] acepack_1.4.1           ellipsis_0.3.2          pkgconfig_2.0.2         XML_3.98-1.9           
# [65] nnet_7.3-12             dbplyr_2.1.1            deldir_0.1-14           locfit_1.5-9.1         
# [69] utf8_1.1.4              tidyselect_1.1.1        rlang_0.4.11            reshape2_1.4.3         
# [73] later_0.8.0             AnnotationDbi_1.38.2    munsell_0.5.0           cellranger_1.1.0       
# [77] tools_3.4.2             cli_3.0.0               generics_0.1.0          RSQLite_2.0            
# [81] broom_0.7.8             yaml_2.2.1              knitr_1.33              bit64_0.9-7            
# [85] fs_1.5.0                robustbase_0.92-8       nlme_3.1-131            mime_0.5               
# [89] xml2_1.3.2              compiler_3.4.2          rstudioapi_0.13         reprex_2.0.0           
# [93] geneplotter_1.54.0      pcaPP_1.9-73            stringi_1.3.1           Matrix_1.2-11          
# [97] vctrs_0.3.8             pillar_1.6.1            LearnBayes_2.15         lifecycle_1.0.0        
# [101] data.table_1.14.0       bitops_1.0-6            httpuv_1.5.0            R6_2.4.0               
# [105] latticeExtra_0.6-28     promises_1.0.1          gridExtra_2.3           codetools_0.2-15       
# [109] lambda.r_1.2            boot_1.3-20             MASS_7.3-47             gtools_3.5.0           
# [113] assertthat_0.2.0        withr_2.4.2             GenomeInfoDbData_0.99.0 mgcv_1.8-22            
# [117] expm_0.999-2            hms_1.1.0               rpart_4.1-11            coda_0.19-1            
# [121] shiny_1.2.0             lubridate_1.7.10        base64enc_0.1-3   
#read in counts
countData <- read.table("Data/GeneExpression_data/host_rtssid_counts_RTonly_BM3.txt")
head(countData)
length(countData[,1])
 #after BM3 9,819
names(countData)=sub(".fastq.trim.sam.counts","",names(countData))
names(countData)

totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts
 # Ssid284  Ssid401  Ssid402  Ssid405  Ssid409  Ssid411  Ssid412  Ssid413  Ssid415  Ssid416  Ssid417 
  # 519518   323215   600920   550094   743124   547549   397676   434612   485601   523617   414490 
 # Ssid418  Ssid429  Ssid430  Ssid435  Ssid436  Ssid438  Ssid440  Ssid444  Ssid445  Ssid446  Ssid447 
  # 613307   496887   534664   564075   430586   651026   331336   503790   478136   613227   321495 
 # Ssid450  Ssid451  Ssid454  Ssid455 Ssid456F  Ssid457  Ssid458  Ssid460  Ssid461  Ssid463  Ssid464 
  # 394210   545799   506455   522783   701497   604132   463545   534075   555605   436165   492648 
 # Ssid467  Ssid470  Ssid471  Ssid472  Ssid476  Ssid477  Ssid495  Ssid497  Ssid498  Ssid499  SsidUNK 
  # 427439   485561   516245   409621   390309   497179   484996   434186   301870   483145   424235
  
min(totalCounts) #301,870
max(totalCounts)  #743,124
mean(totalCounts)

rt <- read.csv("Data/GeneExpression_data/RT_metadata_RTonly.csv")
head(rt)
conditions=data.frame(rt$source, rt$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")

dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~ source+trans) #can only test for the main effects source location and transplant location

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

res<- results(dds)
#############################
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
rld_t=t(rld)
colnames(rld_t)
head(rld_t)
pca <- prcomp(rld_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$trans=conditions$trans
pca_s$source=conditions$source
head(pca_s)

cbPalette <- c("#009E73", "#0072B2", "#D55E00")
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = trans, pch=source, group=trans)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()+
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

adonis(pca_s[,1:2] ~ trans+source, data = pca_s, method='eu', na.rm = TRUE)
#Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# trans      2     402.5 201.274 1.46633 0.06819  0.236
# source     2     147.6  73.808 0.53771 0.02501  0.721
# Residuals 39    5353.3 137.263         0.90681       
# Total     43    5903.4                 1.00000 

#Also try the imbedded plotPCA method within DESeq2
pcaData <- DESeq2::plotPCA(rlog, intgroup = c("source", "trans"), returnData = TRUE) # this one will only pull PC 1&2
cbPalette <- c("#009E73", "#0072B2", "#D55E00")
#pdf("PCA_Host_allgenes_rlog_plotPCA.pdf",height=5,width=6)
ggplot(pcaData, aes(PC1, PC2, color = trans, pch=source, group=trans)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  #xlab(paste0("PC1: ",pc1v,"% variance")) +
  #ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

##try vst normalization
vst_norm=vst(dds, blind=TRUE) 
v=assay(vst_norm)
v_t=t(v)
colnames(v_t)
head(v_t)
pca <- prcomp(v_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$trans=conditions$trans
pca_s$source=conditions$source
head(pca_s)

cbPalette <- c("#009E73", "#0072B2", "#D55E00")
#pdf("PCA_Host_allgenes_vst.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = trans, pch=source, group=trans)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()+
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

pcaData <- DESeq2::plotPCA(vst_norm, intgroup = c("source", "trans"), returnData = TRUE) # this one will only pull PC 1&2
cbPalette <- c("#009E73", "#0072B2", "#D55E00")
#pdf("PCA_Host_allgenes_vst_plotPCA.pdf",height=5,width=6)
ggplot(pcaData, aes(PC1, PC2, color = trans, pch=source, group=trans)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()
# geom_density2d(alpha=.5)+
# geom_polygon(alpha=.2)+
#xlab(paste0("PC1: ",pc1v,"% variance")) +
#ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

##looks the same as rlog- stick with rlog

####################Source NS vs FR
conditions$source<-factor(conditions$source, levels=c("NS","FR"))

resSource1 <- results(dds, contrast=c("source","NS","FR"))
##second term is the "control"

head(resSource1)
#how many FDR < 10%?
table(resSource1$padj<0.1)
# 0.1=5
# 0.05=0
# 0.01=0
head(resSource1)
summary(resSource1)
# out of 9819 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 2, 0.02% 
# LFC < 0 (down)   : 3, 0.031% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 0, 0% 
# (mean count < 3)

nrow(resSource1[resSource1$padj<0.1 & resSource1$log2FoldChange > 0 & !is.na(resSource1$padj),]) 
nrow(resSource1[resSource1$padj<0.1 & resSource1$log2FoldChange <0 & !is.na(resSource1$padj),]) 
#UP in NS 2
#DOWN in NS 3

write.table(resSource1, file="Data/GeneExpression_data/RTsource_NSvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource1=data.frame(resSource1)
head(resSource1)
go_sourceNSvsFR = resSource1 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceNSvsFR)
colnames(go_sourceNSvsFR) <- c("gene", "pval")
head(go_sourceNSvsFR)
write.csv(go_sourceNSvsFR, file="sourceNSvsFR_GO.csv", quote=F, row.names=FALSE)

#########################################################################################################
####################Source NS vs BR
conditions$source<-factor(conditions$source, levels=c("NS","BR"))

resSource2 <- results(dds, contrast=c("source","NS","BR"))
##second term is the "control"
head(resSource2)
#how many FDR < 10%?
table(resSource2$padj<0.1)
# 0.1=0
# 0.05=0
# 0.01=0
head(resSource2)
summary(resSource2)
# out of 9819 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 0, 0% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 0, 0% 
# (mean count < 3)

nrow(resSource2[resSource2$padj<0.1 & resSource2$log2FoldChange > 0 & !is.na(resSource2$padj),]) 
nrow(resSource2[resSource2$padj<0.1 & resSource2$log2FoldChange <0 & !is.na(resSource2$padj),]) 
#UP in NS 0
#DOWN in NS 0

write.table(resSource2, file="Data/GeneExpression_data/RTsource_NSvsBR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource2=data.frame(resSource2)
head(resSource2)
go_sourceNSvsBR = resSource2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceNSvsBR)
colnames(go_sourceNSvsBR) <- c("gene", "pval")
head(go_sourceNSvsBR)
write.csv(go_sourceNSvsBR, file="sourceNSvsBR_GO.csv", quote=F, row.names=FALSE)
#########################################################################################################################
####################Source BR vs FR
conditions$source<-factor(conditions$source, levels=c("BR","FR"))

resSource3 <- results(dds, contrast=c("source","BR","FR"))
##second term is the "control"
head(resSource3)
#how many FDR < 10%?
table(resSource3$padj<0.1)
# 0.1=0
# 0.05=0
# 0.01=0
head(resSource3)
summary(resSource3)
# out of 9819 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 0, 0% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 0, 0% 
# (mean count < 3)

nrow(resSource3[resSource3$padj<0.1 & resSource3$log2FoldChange > 0 & !is.na(resSource3$padj),]) 
nrow(resSource3[resSource3$padj<0.1 & resSource3$log2FoldChange <0 & !is.na(resSource3$padj),]) 
#UP in BR 0
#DOWN in BR 0

write.table(resSource3, file="Data/GeneExpression_data/RTsource_BRvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource3=data.frame(resSource3)
head(resSource3)
go_sourceBRvsFR = resSource3 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceBRvsFR)
colnames(go_sourceBRvsFR) <- c("gene", "pval")
head(go_sourceBRvsFR)
write.csv(go_sourceBRvsFR, file="Data/GeneExpression_data/sourceBRvsFR_GO.csv", quote=F, row.names=FALSE)

####################Transplant NS vs FR
conditions$trans<-factor(conditions$trans, levels=c("NS","FR"))

resTrans1 <- results(dds, contrast=c("trans","NS","FR"))
##second term is the "control"
head(resTrans1)
#how many FDR < 10%?
table(resTrans1$padj<0.01)
# 0.1=78
# 0.05=51
# 0.01=30
head(resTrans1)
summary(resTrans1)
# LFC > 0 (up)     : 20, 0.2% 
# LFC < 0 (down)   : 58, 0.59% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 3994, 41% 
# (mean count < 13)
nrow(resTrans1[resTrans1$padj<0.1 & resTrans1$log2FoldChange > 0 & !is.na(resTrans1$padj),]) 
nrow(resTrans1[resTrans1$padj<0.1 & resTrans1$log2FoldChange <0 & !is.na(resTrans1$padj),]) 
#UP in NS 20
#DOWN in NS 58

write.table(resTrans1, file="Data/GeneExpression_data/RTtransNSvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans1=data.frame(resTrans1)
head(resTrans1)
go_transNSvsFR = resTrans1 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transNSvsFR)
colnames(go_transNSvsFR) <- c("gene", "pval")
head(go_transNSvsFR)
write.csv(go_transNSvsFR, file="Data/GeneExpression_data/transNSvsFR_GO.csv", quote=F, row.names=FALSE)
#########################################################################################################
####################Transplant NS vs BR
conditions$trans<-factor(conditions$trans, levels=c("NS","BR"))

resTrans2=results(dds, contrast=c("trans","NS","BR"))
##second term is the "control"
head(resTrans2)
#how many FDR < 10%?
table(resTrans2$padj<0.11)
# 0.1=103
# 0.05=72
# 0.01=35
head(resTrans2)
summary(resTrans2)
# out of 9819 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 37, 0.38% 
# LFC < 0 (down)   : 64, 0.65% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 1901, 19% 
# (mean count < 6)

nrow(resTrans2[resTrans2$padj<0.1 & resTrans2$log2FoldChange > 0 & !is.na(resTrans2$padj),]) 
nrow(resTrans2[resTrans2$padj<0.1 & resTrans2$log2FoldChange <0 & !is.na(resTrans2$padj),]) 
#UP in NS 37
#DOWN in NS 64

write.table(resTrans2, file="Data/GeneExpression_data/RTtransNSvsBR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans2=data.frame(resTrans2)
head(resTrans2)
go_transNSvsBR = resTrans2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transNSvsBR)
colnames(go_transNSvsBR) <- c("gene", "pval")
head(go_transNSvsBR)
write.csv(go_transNSvsBR, file="Data/GeneExpression_data/transNSvsBR_GO.csv", quote=F, row.names=FALSE)
#########################################################################################################################
####################Transplant BR vs FR
conditions$trans<-factor(conditions$trans, levels=c("BR","FR"))

resTrans3 <- results(dds, contrast=c("trans","BR","FR"))
##second term is the "control"
head(resTrans3)
#how many FDR < 10%?
table(resTrans3$padj<0.01)
# 0.1=283
# 0.05=175
# 0.01=97
head(resTrans3)
summary(resTrans3)
#out of 9819 with nonzero total read count
#adjusted p-value < 0.1
# LFC > 0 (up)     : 115, 1.2% 
# LFC < 0 (down)   : 168, 1.7% 
# outliers [1]     : 8, 0.081% 
# low counts [2]   : 951, 9.7% 
# (mean count < 4)

nrow(resTrans3[resTrans3$padj<0.1 & resTrans3$log2FoldChange > 0 & !is.na(resTrans3$padj),]) 
nrow(resTrans3[resTrans3$padj<0.1 & resTrans3$log2FoldChange <0 & !is.na(resTrans3$padj),]) 
#UP in BR 115
#DOWN in BR 168

write.table(resTrans3, file="Data/GeneExpression_data/RTtransBRvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans3=data.frame(resTrans3)
head(resTrans3)
go_transBRvsFR = resTrans3 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transBRvsFR)
colnames(go_transBRvsFR) <- c("gene", "pval")
head(go_transBRvsFR)
write.csv(go_transBRvsFR, file="Data/GeneExpression_data/transBRvsFR_GO.csv", quote=F, row.names=FALSE)
##############################################################################
#--------------get pvals
head(resSource1)
valSourceNSFR=cbind(resSource1$pvalue, resSource1$padj)
head(valSourceNSFR)
colnames(valSourceNSFR)=c("pval.sNSFR", "padj.sNSFR")
length(valSourceNSFR[,1])
table(complete.cases(valSourceNSFR))

head(resSource2)
valSourceNSBR=cbind(resSource2$pvalue, resSource2$padj)
head(valSourceNSBR)
colnames(valSourceNSBR)=c("pval.sNSBR", "padj.sNSBR")
length(valSourceNSBR[,1])
table(complete.cases(valSourceNSBR))

head(resSource3)
valSourceBRFR=cbind(resSource3$pvalue, resSource3$padj)
head(valSourceBRFR)
colnames(valSourceBRFR)=c("pval.sBRFR", "padj.sBRFR")
length(valSourceBRFR[,1])
table(complete.cases(valSourceBRFR))
 
head(resTrans1)
valTransNSFR=cbind(resTrans1$pvalue, resTrans1$padj)
head(valTransNSFR)
colnames(valTransNSFR)=c("pval.tNSFR", "padj.tNSFR")
length(valTransNSFR[,1])
table(complete.cases(valTransNSFR))

head(resTrans2)
valTransNSBR=cbind(resTrans2$pvalue, resTrans2$padj)
head(valTransNSBR)
colnames(valTransNSBR)=c("pval.tNSBR", "padj.tNSBR")
length(valTransNSBR[,1])
table(complete.cases(valTransNSBR))

head(resTrans3)
valTransBRFR=cbind(resTrans3$pvalue, resTrans3$padj)
head(valTransBRFR)
colnames(valTransBRFR)=c("pval.tBRFR", "padj.tBRFR")
length(valTransBRFR[,1])
table(complete.cases(valTransBRFR))

######-------------make rlogdata (again, used above for PCA) and pvals table
rld=assay(rlog)
head(rld)
conditions=data.frame(rt$source, rt$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")
colnames(rld)=paste(conditions$source, conditions$trans, sep="_")
head(rld)
length(rld[,1]) #9819

rldpvals=cbind(rld,valSourceNSFR, valSourceNSBR, valSourceBRFR, valTransNSFR, valTransNSBR, valTransBRFR)
head(rldpvals)
dim(rldpvals)
# [1]  9819    56
table(complete.cases(rldpvals))
# FALSE  TRUE 
 # 4002  5817

write.csv(rldpvals, "Data/GeneExpression_data/2020_RLDandPVALS_RT_host.csv", quote=F)

#################################################################################
###########################heat map of sample distances for trans
rldpvals <- read.csv(file="Data/GeneExpression_data/2020_RLDandPVALS_RT_host.csv", row.names=1)
gg=read.table("Data/GeneExpression_data/sid_cleaned_iso2gene.tab",sep="\t")

#####trans NS versus BR
head(rldpvals)
colnames(rldpvals)
head(gg)
rld_data=as.data.frame(rldpvals[,c(3,6,7, 9, 10, 11,12,15,17,18,21,23,24,25,29,30,31,33,34,36, 37,39,41,42,43,44,53,54)])
head(rld_data)
colnames(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.tNSBR<=p.val & !is.na(rld_data$padj.tNSBR),]
length(conds[,1])
#101
exp=conds[,1:26]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library("dplyr")
df_all_iso <- explc %>%
rownames_to_column("V1") %>%
left_join(gg) %>%
mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:26]
head(unanno)

df_only_anno <- df_all_iso %>%
filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:26]
head(anno)

df_unanno <- df_all_iso %>%
filter(is.na(V2))

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(anno)

my_sample_col <- data.frame(Source = c("FR", "FR", "FR", "FR", "BR", "NS", "BR", "NS", "BR", "BR", "BR", "FR", "NS", "BR", "BR", "FR", "FR", "NS", "FR", "NS", "BR", "NS", "NS", "NS", "BR"), Transplant = c("BR", "NS", "BR", "BR", "NS", "NS", "BR", "NS", "BR", "NS", "BR","BR","BR","BR","NS", "BR","BR","BR","NS", "BR","BR","BR","NS", "BR", "BR"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  Source = c(NS = "#D55E00", BR = "#009E73", FR = "#0072B2"),
  Transplant = c(NS = "#D55E00", BR = "#009E73", FR = "#0072B2"))

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS4a_NS_BR_DEGS_annotated.pdf",height=8,width=10, onefile = F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames = F, border_color = "NA")
dev.off()

# big heat map of all DEGs, no annotation too
#pdf("Figures/Supplemental_Figures/For_edits/heatmap_NSBR_trans.pdf",height=6,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()
####################################################################################################
#####trans NS versus FR
head(rldpvals)
colnames(rldpvals)
rld_data=as.data.frame(rldpvals[,c(1,2,4,5,6,8,10,11,13,14,15,16,18,19,20,22,26,27,28,29,32,34,35,38,40,41,44,51,52)])
head(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.tNSFR<=p.val & !is.na(rld_data$padj.tNSFR),]
length(conds[,1])
#78
exp=conds[,1:27]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
rownames_to_column("V1") %>%
left_join(gg) %>%
mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:28]
head(unanno)

df_only_anno <- df_all_iso %>%
filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:28]
head(anno)

df_unanno <- df_all_iso %>%
filter(is.na(V2))

colnames(exp)
my_sample_col <- data.frame(Source = c("FR", "NS", "NS","BR", "FR","NS", "BR", "NS", "NS", "NS", "NS", "FR", "BR", "FR", "BR", "NS",  "FR", "BR", "BR", "BR",  "FR", "FR","BR","BR" , "FR", "NS", "NS"), Transplant = c("FR","FR","FR","FR", "NS", "FR","NS", "NS", "FR","FR", "NS","FR", "NS", "FR","FR","FR","FR","FR","FR","NS", "FR", "NS", "FR","FR","FR","NS", "NS"))
row.names(my_sample_col)= colnames(anno)

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS4b_NS_FR_DEGS_annotated.pdf",height=8,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames = F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
#pdf("heatmap_NSFR_trans.pdf",height=6,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()

####################################################################################################
#####trans BR versus FR
head(rldpvals)
colnames(rldpvals)
rld_data=rldpvals[,c(1, 2, 3, 4,5, 7, 8,9,12,13,14,16,17,19, 20, 21,22, 23,24, 25, 26, 27,28, 30, 31, 32, 33, 35,36,37,38,39,40,42,43, 55, 56)]
head(rld_data)
colnames(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.tBRFR<=p.val & !is.na(rld_data$padj.tBRFR),]
length(conds[,1])
#283
exp=conds[,1:35]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
rownames_to_column("V1") %>%
left_join(gg) %>%
mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:36]
head(unanno)

df_only_anno <- df_all_iso %>%
filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:36]
head(anno)

df_unanno <- df_all_iso %>%
filter(is.na(V2))

colnames(exp)
my_sample_col <- data.frame(Source = c("FR","NS","FR","NS","BR","FR","NS","FR","BR", "NS", "NS", "FR", "BR", "FR", "BR", "BR", "NS", "FR", "NS", "BR", "FR", "BR","BR" , "FR","FR","FR", "NS","BR",  "NS", "BR","BR", "NS", "FR","NS", "BR"), Transplant = c("FR","FR","BR","FR","FR", "BR", "FR","BR","BR",  "FR","FR", "FR","BR", "FR","FR","BR","FR","BR","BR","BR","FR","FR","FR","BR","BR","FR","BR","FR","BR","BR","FR","BR","FR","BR", "BR"))
row.names(my_sample_col)= colnames(anno)

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS4c_BR_FR_DEGS_annotated.pdf",height=20,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames = F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
#pdf("heatmap_BRFR_trans.pdf",height=6,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()
################################################################################################################
#####look at the DEG for source NS versus FR
head(rldpvals)
colnames(rldpvals)
rld_data=rldpvals[,c(1, 2, 3,4,6, 7,8,9,11,13,14,15,16,19,22,23,24,26,30,31,32,33,34,36,39,40,41,42,44,45, 46)]
head(rld_data)
colnames(rld_data)
ncol(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.sNSFR<=p.val & !is.na(rld_data$padj.sNSFR),]
length(conds[,1])
#5
exp=conds[,1:29]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
rownames_to_column("V1") %>%
left_join(gg) %>%
mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:30]
head(unanno)

df_only_anno <- df_all_iso %>%
filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:30]
head(anno)
nrow(anno)
anno
df_unanno <- df_all_iso %>%
filter(is.na(V2))

colnames(exp)
my_sample_col <- data.frame(Source = c("FR","NS", "FR", "NS", "FR","FR","NS","FR", "NS", "NS","NS","NS", "FR", "FR", "NS", "FR","NS", "FR", "FR","FR" ,"FR", "NS", "FR", "NS", "NS", "FR", "NS", "NS", "NS"), Transplant = c("FR","FR","BR","FR","NS", "BR","FR","BR", "NS","FR","FR","NS","FR","FR","FR","BR","BR","FR","BR","BR","FR","BR","NS", "BR","BR","FR","NS","BR", "NS"))
row.names(my_sample_col)= colnames(anno)

# cannot make the big heat map of all annotated genes because there are too few to cluster
#Multiple PDZ domain protein is only annotated DEG, whihc was up in NS relative to FR 

# big heat map of all DEGs, no annotation too

####################################################################################################
####################################################################################################
#####source NS versus BR
head(rldpvals)
conds=rldpvals[rldpvals$padj.sNSBR<=p.val & !is.na(rldpvals$padj.sNSBR),]
length(conds[,1])
#no DEGS

####################################################################################################
#####source BR versus FR
head(rldpvals)
conds=rldpvals[rldpvals$padj.sBRFR<=p.val & !is.na(rldpvals$padj.sBRFR),]
length(conds[,1])
#0 DEGS

###########################################################################################
# VENN Diagram to include both up and down regulated genes in common for source
# install.packages('VennDiagram', dependencies=TRUE, repos='http://cran.us.r-project.org')
library(VennDiagram)
head(resSource1)
head(resSource2)
head(resSource3)

pSource1up=row.names(resSource1[resSource1$padj<0.1 & !is.na(resSource1$padj) & resSource1$log2FoldChange>0,])
length(pSource1up) #2
pSource1down=row.names(resSource2[resSource1$padj<0.1 & !is.na(resSource1$padj) & resSource1$log2FoldChange<0,])
length(pSource1down) #3
pSource2up=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj) & resSource2$log2FoldChange>0,])
length(pSource2up) #0
pSource2down=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj) & resSource2$log2FoldChange<0,])
length(pSource2down) #0
pSource3up=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj) & resSource3$log2FoldChange>0,])
length(pSource3up) #0
pSource3down=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj) & resSource3$log2FoldChange<0,])
length(pSource3down) #0

pSource1=row.names(resSource1[resSource1$padj<0.1 & !is.na(resSource1$padj),])
pSource2=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj),])
pSource3=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj),])

#UP
pdegs1_up=union(pSource1up,pSource2up)
pdegs1_up=union(pdegs1_up,pSource3up)
length(pdegs1_up)
#2

#DOWN
pdegs1_down=union(pSource1down,pSource2down)
pdegs1_down=union(pdegs1_down,pSource3down)
length(pdegs1_down)
#3

#ALL
pdegs1=union(pSource1,pSource2)
pdegs1=union(pdegs1,pSource3)
length(pdegs1)
#5

###do UP, DOWN, ALL
candidates=list("NSvsFRdown"=pSource1down, "NSvsBRdown"=pSource2down, "BRvsFRdown"=pSource3down)
#quartz()
#prettyvenn=venn.diagram(
#  x = candidates,
#  filename=NULL,
#  col = "transparent",
#  fill = c("coral2", "forestgreen", "royalblue1"),
#  alpha = 0.5,
#  label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
#  cex = 2.5,
#  fontfamily = "sans",
#  fontface = "bold",
#  cat.default.pos = "text",
#  cat.col = c("darkred", "darkgreen", "blue4"),
#  cat.cex = 2.5,
#  cat.fontfamily = "sans",
#  cat.dist = c(0.08, 0.08, 0.03),
#  cat.pos = 1
#);
#grid.draw(prettyvenn)

# VENN Diagram to include both up and down regulated genes in common for transplant
head(resTrans1)
head(resTrans2)
head(resTrans3)

pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #20
pTrans1down=row.names(resTrans2[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #58
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #37
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #64
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #115
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #168

pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])

#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
length(pdegs1_up)
#160

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
length(pdegs1_down)
#245

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
length(pdegs1)
#364

###do UP, DOWN, ALL
candidates=list("NSvsFRup"=pTrans1up, "NSvsBRup"=pTrans2up, "BRvsFRup"=pTrans3up)
#quartz()
#prettyvenn=venn.diagram(
#  x = candidates,
#  filename=NULL,
#  col = "transparent",
#  fill = c("coral2", "forestgreen", "royalblue1"),
#  alpha = 0.5,
#  label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
#  cex = 2.5,
#  fontfamily = "sans",
#  fontface = "bold",
#  cat.default.pos = "text",
#  cat.col = c("darkred", "darkgreen", "blue4"),
#  cat.cex = 2.5,
#  cat.fontfamily = "sans",
#  cat.dist = c(0.08, 0.08, 0.03),
#  cat.pos = 1
#);
#grid.draw(prettyvenn)


# CANONICAL CORRESPONDENCE ANALYSIS OF GENE EXPRESSION based on raw counts
##packages used in plasticity analysis. need to restart R to avoid DLL issue
library("ggpubr")
library("tidyverse")
library("dplyr")
library("performance")
library("finalfit")
library("see")
library("vegan")
library("adegenet") 
library("WGCNA") # for labels2colors

sessionInfo()
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] WGCNA_1.61            fastcluster_1.1.24    dynamicTreeCut_1.63-1 adegenet_2.1.0       
# [5] ade4_1.7-8            vegan_2.5-4           lattice_0.20-35       permute_0.9-4        
# [9] see_0.7.0             finalfit_1.0.4        performance_0.9.0     forcats_0.5.1        
# [13] stringr_1.4.0         dplyr_1.0.7           purrr_0.3.4           readr_1.4.0          
# [17] tidyr_1.1.3           tibble_3.1.2          tidyverse_1.3.1       ggpubr_0.2           
# [21] magrittr_2.0.1        ggplot2_3.3.5        
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.1.2       Hmisc_4.1-0           plyr_1.8.4           
# [5] igraph_1.1.2          sp_1.2-7              splines_3.4.2         robust_0.4-18        
# [9] digest_0.6.18         foreach_1.4.3         htmltools_0.3.6       GO.db_3.4.1          
# [13] gdata_2.18.0          fansi_0.4.0           checkmate_1.8.5       memoise_1.1.0        
# [17] fit.models_0.5-14     cluster_2.0.6         doParallel_1.0.11     modelr_0.1.8         
# [21] matrixStats_0.52.2    gmodels_2.16.2        colorspace_1.4-0      rrcov_1.4-3          
# [25] blob_1.2.1            rvest_1.0.0           xfun_0.24             haven_2.4.1          
# [29] pan_1.6               crayon_1.4.1          jsonlite_1.7.2        lme4_1.1-21          
# [33] impute_1.50.1         survival_2.41-3       iterators_1.0.8       ape_5.5              
# [37] glue_1.4.2            gtable_0.2.0          seqinr_3.4-5          DEoptimR_1.0-8       
# [41] BiocGenerics_0.22.1   jomo_2.6-7            scales_1.0.0          mvtnorm_1.0-7        
# [45] DBI_1.1.1             Rcpp_1.0.0            htmlTable_1.11.1      xtable_1.8-2         
# [49] spData_0.2.6.9        foreign_0.8-69        bit_1.1-12            spdep_0.7-4          
# [53] preprocessCore_1.38.1 Formula_1.2-2         stats4_3.4.2          htmlwidgets_1.3      
# [57] httr_1.4.2            RColorBrewer_1.1-2    acepack_1.4.1         ellipsis_0.3.2       
# [61] mice_3.4.0            pkgconfig_2.0.2       nnet_7.3-12           dbplyr_2.1.1         
# [65] deldir_0.1-14         utf8_1.1.4            tidyselect_1.1.1      rlang_0.4.11         
# [69] reshape2_1.4.3        later_0.8.0           AnnotationDbi_1.38.2  munsell_0.5.0        
# [73] cellranger_1.1.0      tools_3.4.2           cli_3.0.0             generics_0.1.0       
# [77] RSQLite_2.0           broom_0.7.8           yaml_2.2.1            knitr_1.33           
# [81] bit64_0.9-7           fs_1.5.0              robustbase_0.92-8     mitml_0.4-3          
# [85] nlme_3.1-131          mime_0.5              xml2_1.3.2            compiler_3.4.2       
# [89] rstudioapi_0.13       reprex_2.0.0          pcaPP_1.9-73          stringi_1.3.1        
# [93] Matrix_1.2-11         nloptr_1.2.1          vctrs_0.3.8           pillar_1.6.1         
# [97] LearnBayes_2.15       lifecycle_1.0.0       data.table_1.14.0     insight_0.17.1       
# [101] httpuv_1.5.0          R6_2.4.0              latticeExtra_0.6-28   promises_1.0.1       
# [105] gridExtra_2.3         IRanges_2.10.5        codetools_0.2-15      boot_1.3-20          
# [109] MASS_7.3-47           gtools_3.5.0          assertthat_0.2.0      withr_2.4.2          
# [113] S4Vectors_0.14.7      mgcv_1.8-22           expm_0.999-2          parallel_3.4.2       
# [117] hms_1.1.0             grid_3.4.2            rpart_4.1-11          coda_0.19-1          
# [121] minqa_1.2.4           Biobase_2.36.2        shiny_1.2.0           lubridate_1.7.10     
# [125] base64enc_0.1-3  

cc <- read.table("Data/GeneExpression_data/host_rtssid_counts_RTonly_BM3.txt")
head(cc)
names(cc)=sub(".fastq.trim.sam.counts","",names(cc))
names(cc)
ncol(cc)

mns=apply(cc,1,mean)
table(mns>1)
cc=cc[mns>1,]
cct=t(cc) # in vegan, species (genes) must be columns, and sites (samples) must be rows
head(cct)
names(cc)
length(cc)

rt=read.csv("Data/GeneExpression_data/RT_metadata_RTonly.csv")
conditions=data.frame(rt$source, rt$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")

# normalizing to total counts
cct=decostand(cct,method="total")
# log-transform leaving zeroes as zeroes
cct=decostand(cct,method="log")

# canonical correspondence analysis (chi-square distances) - prepare to wait a while
ccaa=cca(cct~., conditions)
head(ccaa)
summary(ccaa)
#Accumulated constrained eigenvalues
# Importance of components:
#   CCA1      CCA2      CCA3      CCA4
# Eigenvalue            0.0008343 0.0007439 0.0006873 0.0006044
# Proportion Explained  0.2907176 0.2592047 0.2394874 0.2105904
# Cumulative Proportion 0.2907176 0.5499222 0.7894096 1.0000000

pt_colors <- ifelse(conditions$source == 'FR', '#0072B2',
                    ifelse(conditions$source == 'NS', '#D55E00', '#009E73'))
pt_transcol <- ifelse(conditions$trans == 'FR', '#0072B2',
                    ifelse(conditions$trans == 'NS', '#D55E00', '#009E73'))
shape_source <- ifelse(conditions$source == 'FR', 15,
                    ifelse(conditions$source == 'NS', 16, 17))
shape_trans <- ifelse(conditions$trans == 'FR', 15,
                       ifelse(conditions$trans == 'NS', 16, 17))

pdf("Figures/Supplemental_Figures/For_edits/Fig3b_CCA_Host_trans_final.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites", type="n", ylab="CCA2 (25.9%)", xlab="CCA1 (29.1%)", xlim=range(-2.5,3), ylim=range(-2.5,2.1))
points(ccaa,choices=c(1,2),col=pt_colors, pch=shape_trans)
ordispider(ccaa,groups=conditions$trans,col="grey50")
ordiellipse(ccaa, groups=conditions$trans, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 1.5, cex = 0.8, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("grey","grey","grey"), pch = c(16,17,15),
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.13, 0.005), title="Transplant")
dev.off()

pdf("Figures/Supplemental_Figures/For_edits/Fig3a_CCA_host_source_final.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites",type="n",  ylab="CCA2 (25.9%)", xlab="CCA1 (29.1%)", xlim=range(-2.5,3), ylim=range(-2.5,2.1))
points(ccaa,choices=c(1,2),col=pt_colors,pch=shape_trans)
ordispider(ccaa,groups=conditions$source,col="grey50")
ordiellipse(ccaa, groups=conditions$source, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 1.5, cex = 0.8, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("grey","grey","grey"), pch = c(16,17,15),
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.13, 0.005), title="Transplant")
dev.off()

#----------------------------
# ANOVA
anova(ccaa,by = "term") 
# Permutation test for cca under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# Model: cca(formula = cct ~ source + trans, data = conditions)
# Df ChiSquare      F Pr(>F)  
# source    2 0.0013555 1.0331  0.229  
# trans     2 0.0015144 1.1542  0.015 *
#   Residual 39 0.0255842 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## -------------------------------------------------------- Colleen's additions below for plasticity:
set.seed(16)
#### I decided it will be cleaner if you source these functions from another script so you will need to do that now! Just save it to your working directory or add the path to it in the source function:
source("Code/CCAplast_function.R")
#quartz()
screeplot(ccaa)
## create dataframe of CCAs with experimental conditions
ccaa_full <- data.frame(conditions,
                        trt = paste(conditions$source, conditions$trans, sep = "_"),
                        scores(ccaa, choices = c(1:4))$sites)
# create dataframes for each source location
ccaa_fr <- ccaa_full %>% filter(trans == "FR")
ccaa_br <- ccaa_full %>% filter(trans == "BR")
ccaa_ns <- ccaa_full %>% filter(trans == "NS")

#### Calculate plasticity for each source location with custom function:
## Right now, these are calculating the plasticity based on all CCAs
# If you want to change that, you can specify it by unhashing 'num_cca' and specifying how many CCAs to include

plast_FR <-  CCAplast(cca = ccaa_fr[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_fr[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "FR_FR", # control level of the treatment. If blank, a control mean per control level is assumed
                      keep = "yes") # change to 'yes' if you want to retain the control level values

plast_BR <-  CCAplast(cca = ccaa_br[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_br[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "BR_BR", # control level of the treatment. If blank, a control mean per control level is assumed   
                      keep = "yes") # change to 'yes' if you want to retain the control level values

plast_NS <-  CCAplast(cca = ccaa_ns[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_ns[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "NS_NS", # control level of the treatment. If blank, a control mean per control level is assumed
                      keep = "yes") # change to 'yes' if you want to retain the control level values

## combine all three plasticity dataframes into one
ccaa_plast <- rbind(plast_FR, plast_BR, plast_NS)


##############
## Selecting best-fit model of GLM via AIC
# need to first filter out the single FR -> NS sample

## run a full model
plast_mod <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_plast)
check_model(plast_mod)
## data look fairly normal so proceeding with GLM here ##
## fit different model combos
plast_mod <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_plast) # best-fit model
plast_mod2 <- glm(dist ~ source + trans, family = Gamma(link = "log"), data = ccaa_plast)
plast_mod3 <- glm(dist ~ trans, family = Gamma(link = "log"), data = ccaa_plast)
plast_mod4 <- glm(dist ~ source, family = Gamma(link = "log"), data = ccaa_plast)
plast_mod5 <- glm(dist ~ trt, family = Gamma(link = "log"), data = ccaa_plast)
## check for best-fit model
compare_performance(plast_mod, plast_mod2, plast_mod3, plast_mod4, plast_mod5)
# Name       | Model |     AIC | AIC weights |     BIC | BIC weights | Nagelkerke's R2 |  RMSE | Sigma
# ----------------------------------------------------------------------------------------------------
# plast_mod  |   glm |  46.845 |       0.500 |  64.687 |       0.500 |           0.833 | 0.316 | 0.478
# plast_mod2 |   glm | 102.088 |     < 0.001 | 112.794 |     < 0.001 |           0.146 | 0.804 | 0.891
# plast_mod3 |   glm |  99.632 |     < 0.001 | 106.769 |     < 0.001 |           0.108 | 0.770 | 0.883
# plast_mod4 |   glm | 102.398 |     < 0.001 | 109.535 |     < 0.001 |           0.034 | 0.800 | 0.908
# plast_mod5 |   glm |  46.845 |       0.500 |  64.687 |       0.500 |           0.833 | 0.316 | 0.478

## Best-fit GLMM with Gamma log link
plast_glm <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_plast)
summary(plast_glm) # summary output

# Call:
#   glm(formula = dist ~ source * trans, family = Gamma(link = "log"), 
#       data = ccaa_plast)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.42485  -0.22767  -0.05226   0.23247   0.79365  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -0.6692     0.1731  -3.866 0.000460 ***
#   sourceFR           1.1464     0.2448   4.683 4.17e-05 ***
#   sourceNS           0.6369     0.2568   2.481 0.018072 *  
#   transFR            1.1503     0.2448   4.699 3.97e-05 ***
#   transNS            0.0778     0.2998   0.259 0.796787    
# sourceFR:transFR  -3.0430     0.3462  -8.789 2.22e-10 ***
#   sourceNS:transFR  -0.2569     0.3548  -0.724 0.473857    
# sourceFR:transNS   0.1188     0.4580   0.259 0.796844    
# sourceNS:transNS  -1.5480     0.4133  -3.746 0.000647 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.1798039)
# 
# Null deviance: 34.6034  on 43  degrees of freedom
# Residual deviance:  8.0074  on 35  degrees of freedom
# AIC: 46.845
# 
# Number of Fisher Scoring iterations: 5
### Bootstrapping with custom function
# This is a bootstrap function that will sample WITH replacement (explained step by step within the function)
bootFUN <- function(model, newdata) {
  if(class(model)[1] == "negbin"){
    nr <- nrow(model[["model"]]) # count number of rows in the model (# of observations)
    data <- model[["model"]] # pull data from model
  } else {
    nr <- nrow(model$data) # count number of rows in the model (# of observations)
    data <- model$data # pull data from model
  }
  data2 <- data.frame(data[sample(1:nr, replace = TRUE), ]) # random sample of numbers from 1 - nr (# of rows) to select data from
  update <- update(model, data = data2) # rerun the model using the 'new' dataset from the random smapling above
  predict(update, newdata, type = "response", level = 0.95, allow.new.levels = TRUE) # predicts response variable with model and updated data 
}
bootnum <- 1500 # this is the number of bootstrap replicates. I set it to 1500 here but we should select a number between 999 and 9999

# construct new dataframe and then run bootstrapping
newdata_boot <- data.frame(coral = rownames(ccaa_plast),trt = ccaa_plast$trt, source = ccaa_plast$source, trans = ccaa_plast$trans, dist = ccaa_plast$dist)
boot <- replicate(bootnum, bootFUN(model = plast_glm, newdata = newdata_boot)) # run the bootFUN the number (bootnum) of times specified and save as matrix

# Calculate the mean, 95% lowerCI, and 95% upperCI from the boot matrix and add it to dataframe (but as a new name)
plast_boot <- cbind(ccaa_plast, as.data.frame(t(apply(boot, 1, function(x) c(quantile(x, c(0.025, 0.5, 0.975)))))))
colnames(plast_boot)[5:7] <- c("lowerci", "estimate", "upperci") # rename mean/CI columns


##############
plast_boot$source <- relevel(plast_boot$source, ref = "NS")
plast_boot$trans <- relevel(plast_boot$trans, ref = "NS")

# plot of plasticity (host)

ggplot(data = plast_boot, aes(x = trans, colour = source)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 0, 0))), line = "none") +  
  geom_point(aes(y = dist), position = position_jitterdodge(dodge.width = 0.3), alpha = 0.4, size = 1, shape = 1) +
  geom_linerange(aes(ymin = lowerci, ymax = upperci, colour = source), size = 0.7, position = position_dodge(width = 0.3)) +
  geom_point(aes(y = estimate), position = position_dodge(width = 0.3), size = 2) +
  scale_color_manual("Source", values = c("#D55E00", "#009E73", "#0072B2")) +
  labs(x = "Transplant Location", y = "Plasticity")
ggsave("Figures/Supplemental_Figures/For_edits/Fig3c_CCA_host_plasticity_CI.pdf", height = 3, width = 3, useDingbats = FALSE)

ccaa_plast_sum <- ccaa_plast %>%
  group_by(trans, source) %>%
  summarise(mean = mean(dist, na.rm = TRUE),
            sd = sd(dist))
#relevel
ccaa_plast$source <- relevel(ccaa_plast$source, ref = "NS")
ccaa_plast$trans <- relevel(ccaa_plast$trans, ref = "NS")
# plot of plasticity
ggplot(data = ccaa_plast, aes(x = trans, colour = source)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 0, 0))), line = "none") +  
  geom_point(aes(y = dist), position =position_jitterdodge(dodge.width = 0.3), alpha = 0.4, size = 1, shape = 1) +
  geom_errorbar(data = ccaa_plast_sum, aes(x = trans, ymin = mean - sd, ymax = mean + sd, colour = source), size = 0.7, width = 0.0, position = position_dodge(width = 0.2)) +
  geom_point(data = ccaa_plast_sum, aes(y = mean), position = position_dodge(width = 0.2), size = 2) +
  scale_color_manual("Source", values = c("#D55E00", "#009E73", "#0072B2")) +
  labs(x = "Transplant Location", y = "Plasticity")
#ggsave("Figures/Supplemental_Figures/For_edits/Fig3c_CCA_host_plasticity_SD.pdf", height = 3, width = 3, useDingbats = FALSE)


###Data were then put into GO_MWU for GO enrichment analysis to associate function








##### Ssid reciprocal transplant RNAseq analysis
### algal symbiont data - all assumed to be Cladocopium reads
#setwd("~/Dropbox/UNC/RT_Paper/2022_Final_RT/GE")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("vegan")
library("tidyverse")
library("pheatmap")
library("VennDiagram")
library("adegenet") 
library("WGCNA")

sessionInfo()
#same as above with host data

#read in counts
# countData <- read.table("Data/GeneExpression_data/sym_rtssid_counts.txt")
countData <- read.table("Data/GeneExpression_data/sym_rtssid_counts_BM3.txt")
head(countData)
length(countData[,1])
#23944  #after BM3 4750
names(countData)=sub(".fastq.trim.sam.counts","",names(countData))
names(countData)

totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts
# Ssid284  Ssid401  Ssid402  Ssid405  Ssid409  Ssid411  Ssid412  Ssid413  Ssid415  Ssid416  Ssid417 
# 63016    38890   168102    67775    79255    59781    89074    68982    65403    54596    49852 
# Ssid418  Ssid429  Ssid430  Ssid435  Ssid436  Ssid438  Ssid440  Ssid444  Ssid445  Ssid446  Ssid447 
# 55579    78184    79425    61490    60565    96726    41859    72573    93115    72815    33569 
# Ssid450  Ssid451  Ssid454  Ssid455 Ssid456F  Ssid457  Ssid458  Ssid460  Ssid461  Ssid463  Ssid464 
# 48953    61487    55569    61426    63886    77722    52855    68271    90576    54794    64890 
# Ssid467  Ssid470  Ssid471  Ssid472  Ssid476  Ssid477  Ssid495  Ssid497  Ssid498  Ssid499  SsidUNK 
# 29442    69572    49637   106843    44312    77713    91491    56860    23375    50486    43552

min(totalCounts) #23,375
max(totalCounts)  #168,102

###remove individuals that have very low counts 
rt=read.csv("Data/GeneExpression_data/RT_metadata_RTonly.csv")
conditions=data.frame(rt$source, rt$transplant)
nrow(conditions) 
ncol(countData)
names(conditions)=c("source", "trans")
###4 outliers detected (401,447, 467, 498)
#401: 38890
#447: 33569
#467: 29442
#498: 23375

#remove samples with <40k reads
countData$Ssid401 <- NULL
countData$Ssid447 <- NULL
countData$Ssid467 <- NULL
countData$Ssid498 <- NULL
ncol(countData)

totalCounts=colSums(countData)
totalCounts
min(totalCounts) #41859
max(totalCounts)  #168,102

rt2=rt[!(rt$sample %in% c("Ssid401", "Ssid447", "Ssid467", "Ssid498")), ]
conditions=data.frame(rt2$source, rt2$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")

dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~ source+trans) #can only test for the main effects source location and transplant location
dds_symb <- DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~ source+trans) 

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

res<- results(dds)

#############################
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
rld_t=t(rld)
colnames(rld_t)
head(rld_t)
pca <- prcomp(rld_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$trans=conditions$trans
pca_s$source=conditions$source
head(pca_s)

cbPalette <- c("#009E73", "#0072B2", "#D55E00")
#pdf("PCA_Sym_allgenes_source.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = source, pch = trans, group=source)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()+
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

adonis(pca_s[,1:2] ~ trans+source, data = pca_s, method='eu', na.rm = TRUE)
#Permutation: free
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# trans      2    121.35  60.675  1.5749 0.05992  0.190    
# source     2    555.54 277.769  7.2097 0.27429  0.001 ***
#   Residuals 35   1348.45  38.527         0.66579           
# Total     39   2025.34                 1.00000           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

####################Source NS vs FR
conditions$source<-factor(conditions$source, levels=c("NS","FR"))

resSource1 <- results(dds, contrast=c("source","NS","FR"))
head(resSource1)
#how many FDR < 10%?
table(resSource1$padj<0.05)
# 0.1=29
# 0.05=17
# 0.01=6
head(resSource1)
summary(resSource1)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 16, 0.34%
# LFC < 0 (down)     : 13, 0.27%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 2577, 54%
# (mean count < 6)

nrow(resSource1[resSource1$padj<0.1 & resSource1$log2FoldChange > 0 & !is.na(resSource1$padj),]) 
nrow(resSource1[resSource1$padj<0.1 & resSource1$log2FoldChange <0 & !is.na(resSource1$padj),]) 
#UP in NS 16
#DOWN in NS 13

write.table(resSource1, file="Data/GeneExpression_data/sym_RTsource_NSvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource1=data.frame(resSource1)
head(resSource1)
go_sourceNSvsFR = resSource1 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceNSvsFR)
colnames(go_sourceNSvsFR) <- c("gene", "pval")
head(go_sourceNSvsFR)
write.csv(go_sourceNSvsFR, file="Data/GeneExpression_data/source_sym_NSvsFR_GO.csv", quote=F, row.names=FALSE)
#########################################################################################################
####################Source NS vs BR
conditions$source<-factor(conditions$source, levels=c("NS","BR"))

resSource2 <- results(dds, contrast=c("source","NS","BR"))
head(resSource2)
#how many FDR < 10%?
table(resSource2$padj<0.01)
# 0.1=9
# 0.05=6
# 0.01=3
head(resSource2)
summary(resSource2)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5, 0.11%
# LFC < 0 (down)     : 4, 0.084%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 1195, 25%
# (mean count < 4)

nrow(resSource2[resSource2$padj<0.1 & resSource2$log2FoldChange > 0 & !is.na(resSource2$padj),]) 
nrow(resSource2[resSource2$padj<0.1 & resSource2$log2FoldChange <0 & !is.na(resSource2$padj),]) 
#UP in NS 5
#DOWN in NS 4

write.table(resSource2, file="Data/GeneExpression_data/sym_RTsource_NSvsBR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource2=data.frame(resSource2)
head(resSource2)
go_sourceNSvsBR = resSource2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceNSvsBR)
colnames(go_sourceNSvsBR) <- c("gene", "pval")
head(go_sourceNSvsBR)
write.csv(go_sourceNSvsBR, file="Data/GeneExpression_data/source_sym_NSvsBR_GO.csv", quote=F, row.names=FALSE)

#########################################################################################################################
####################Source BR vs FR
conditions$source<-factor(conditions$source, levels=c("BR","FR"))

resSource3 <- results(dds, contrast=c("source","BR","FR"))
head(resSource3)
#how many FDR < 10%?
table(resSource3$padj<0.1)
# 0.1=73
# 0.05=41
# 0.01=7
head(resSource3)
summary(resSource3)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 28, 0.59%
# LFC < 0 (down)     : 45, 0.95%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 3313, 70%
# (mean count < 9)

nrow(resSource3[resSource3$padj<0.1 & resSource3$log2FoldChange > 0 & !is.na(resSource3$padj),]) 
nrow(resSource3[resSource3$padj<0.1 & resSource3$log2FoldChange <0 & !is.na(resSource3$padj),]) 
#UP in BR 28
#DOWN in BR 45

write.table(resSource3, file="Data/GeneExpression_data/sym_RTsource_BRvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resSource3=data.frame(resSource3)
head(resSource3)
go_sourceBRvsFR = resSource3 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_sourceBRvsFR)
colnames(go_sourceBRvsFR) <- c("gene", "pval")
head(go_sourceBRvsFR)
write.csv(go_sourceBRvsFR, file="Data/GeneExpression_data/source_sym_BRvsFR_GO.csv", quote=F, row.names=FALSE)

####################Transplant NS vs FR
conditions$trans<-factor(conditions$trans, levels=c("NS","FR"))

resTrans1 <- results(dds, contrast=c("trans","NS","FR"))
head(resTrans1)
#how many FDR < 10%?
table(resTrans1$padj<0.1)
# 0.1=0
# 0.05=0
# 0.01=0
head(resTrans1)
summary(resTrans1)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 0, 0%
# (mean count < 1)
nrow(resTrans1[resTrans1$padj<0.1 & resTrans1$log2FoldChange > 0 & !is.na(resTrans1$padj),]) 
nrow(resTrans1[resTrans1$padj<0.1 & resTrans1$log2FoldChange <0 & !is.na(resTrans1$padj),]) 
#UP in NS 0
#DOWN in NS 0

write.table(resTrans1, file="Data/GeneExpression_data/sym_RTtransNSvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans1=data.frame(resTrans1)
head(resTrans1)
go_transNSvsFR = resTrans1 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transNSvsFR)
colnames(go_transNSvsFR) <- c("gene", "pval")
head(go_transNSvsFR)
write.csv(go_transNSvsFR, file="Data/GeneExpression_data/trans_sym_NSvsFR_GO.csv", quote=F, row.names=FALSE)

#########################################################################################################
####################Transplant NS vs BR
conditions$trans<-factor(conditions$trans, levels=c("NS","BR"))

resTrans2=results(dds, contrast=c("trans","NS","BR"))
head(resTrans2)
#how many FDR < 10%?
table(resTrans2$padj<0.1)
# 0.1=0
# 0.05=0
# 0.01=0
head(resTrans2)
summary(resTrans2)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 0, 0%
# (mean count < 1)

nrow(resTrans2[resTrans2$padj<0.1 & resTrans2$log2FoldChange > 0 & !is.na(resTrans2$padj),]) 
nrow(resTrans2[resTrans2$padj<0.1 & resTrans2$log2FoldChange <0 & !is.na(resTrans2$padj),]) 
#UP in NS 0
#DOWN in NS 0

write.table(resTrans2, file="Data/GeneExpression_data/sym_RTtransNSvsBR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans2=data.frame(resTrans2)
head(resTrans2)
go_transNSvsBR = resTrans2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transNSvsBR)
colnames(go_transNSvsBR) <- c("gene", "pval")
head(go_transNSvsBR)
write.csv(go_transNSvsBR, file="Data/GeneExpression_data/trans_sym_NSvsBR_GO.csv", quote=F, row.names=FALSE)
#########################################################################################################################
####################Transplant BR vs FR
conditions$trans<-factor(conditions$trans, levels=c("BR","FR"))

resTrans3 <- results(dds, contrast=c("trans","BR","FR"))
head(resTrans3)
#how many FDR < 10%?
table(resTrans3$padj<0.1)
# 0.1=3
# 0.05=0
# 0.01=0
head(resTrans3)
summary(resTrans3)
# out of 4750 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 3, 0.063%
# outliers [1]       : 2, 0.042%
# low counts [2]     : 0, 0%
# (mean count < 1)

nrow(resTrans3[resTrans3$padj<0.1 & resTrans3$log2FoldChange > 0 & !is.na(resTrans3$padj),]) 
nrow(resTrans3[resTrans3$padj<0.1 & resTrans3$log2FoldChange <0 & !is.na(resTrans3$padj),]) 
#UP in BR 0
#DOWN in BR 3

write.table(resTrans3, file="Data/GeneExpression_data/sym_RTtransBRvsFR.txt", quote=F, sep="\t")

##make the GO table for MWU
resTrans3=data.frame(resTrans3)
head(resTrans3)
go_transBRvsFR = resTrans3 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_transBRvsFR)
colnames(go_transBRvsFR) <- c("gene", "pval")
head(go_transBRvsFR)
write.csv(go_transBRvsFR, file="Data/GeneExpression_data/trans_sym_BRvsFR_GO.csv", quote=F, row.names=FALSE)


##############################################################################
#--------------get pvals
head(resSource1)
valSourceNSFR=cbind(resSource1$pvalue, resSource1$padj)
head(valSourceNSFR)
colnames(valSourceNSFR)=c("pval.sNSFR", "padj.sNSFR")
length(valSourceNSFR[,1])
table(complete.cases(valSourceNSFR))

head(resSource2)
valSourceNSBR=cbind(resSource2$pvalue, resSource2$padj)
head(valSourceNSBR)
colnames(valSourceNSBR)=c("pval.sNSBR", "padj.sNSBR")
length(valSourceNSBR[,1])
table(complete.cases(valSourceNSBR))

head(resSource3)
valSourceBRFR=cbind(resSource3$pvalue, resSource3$padj)
head(valSourceBRFR)
colnames(valSourceBRFR)=c("pval.sBRFR", "padj.sBRFR")
length(valSourceBRFR[,1])
table(complete.cases(valSourceBRFR))

head(resTrans1)
valTransNSFR=cbind(resTrans1$pvalue, resTrans1$padj)
head(valTransNSFR)
colnames(valTransNSFR)=c("pval.tNSFR", "padj.tNSFR")
length(valTransNSFR[,1])
table(complete.cases(valTransNSFR))

head(resTrans2)
valTransNSBR=cbind(resTrans2$pvalue, resTrans2$padj)
head(valTransNSBR)
colnames(valTransNSBR)=c("pval.tNSBR", "padj.tNSBR")
length(valTransNSBR[,1])
table(complete.cases(valTransNSBR))

head(resTrans3)
valTransBRFR=cbind(resTrans3$pvalue, resTrans3$padj)
head(valTransBRFR)
colnames(valTransBRFR)=c("pval.tBRFR", "padj.tBRFR")
length(valTransBRFR[,1])
table(complete.cases(valTransBRFR))

######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
conditions=data.frame(rt2$source, rt2$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")
colnames(rld)=paste(conditions$source, conditions$trans, sep="_")
head(rld)
length(rld[,1]) #4750

rldpvals=cbind(rld,valSourceNSFR, valSourceNSBR, valSourceBRFR, valTransNSFR, valTransNSBR, valTransBRFR)
head(rldpvals)
dim(rldpvals)
# [1]  4750    52
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 3315  1435

write.csv(rldpvals, "Data/GeneExpression_data/2020_RLDandPVALS_RT_sym.csv", quote=F)

#################################################################################
###########################heat map of sample distances for trans
rldpvals <- read.csv(file="Data/GeneExpression_data/2020_RLDandPVALS_RT_sym.csv", row.names=1)
gg=read.table("Data/GeneExpression_data/sym_feb_iso2gene.tab",sep="\t")

#####trans NS versus BR
head(rldpvals)
head(gg)
colnames(rldpvals)
rld_data=as.data.frame(rldpvals[,c(2,5, 6, 8, 9, 10, 11, 14, 16,17, 20,21,22,23,27,28, 29, 31,33,34,36, 38,39, 49, 50)])
head(rld_data)
colnames(rld_data)
p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.tNSBR<=p.val & !is.na(rld_data$padj.tNSBR),]
length(conds[,1])
#0

####source NS versus BR
colnames(rldpvals)
rld_data=as.data.frame(rldpvals[,c(3,4, 7, 9, 10, 11, 12, 13, 14, 16,17,19, 20,22,23,25,26,27, 31,32,33,34,35,36, 38,39,40, 43,44)])
head(rld_data)
colnames(rld_data)
p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.sNSBR<=p.val & !is.na(rld_data$padj.sNSBR),]
length(conds[,1])
#9
head(conds)

exp=conds[,1:27]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library(tidyverse)
df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:28]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:28]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(anno)

my_sample_col <- data.frame(Source = c("NS", "BR", "NS", "BR", "NS", "BR", "NS","NS","NS", "BR", "BR", "BR", "BR", "NS", "BR", "BR", "BR","BR", "NS", "BR", "NS",  "BR","BR", "NS", "NS", "BR", "NS"), 
                            Transplant = c("FR","FR","FR", "NS",  "NS", "BR","FR","FR", "NS", "BR", "NS", "FR","BR","BR","BR","FR","FR","NS", "BR", "FR","BR", "BR", "FR","BR", "NS", "BR", "NS"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  Source = c(NS = "#D55E00", BR = "#009E73", FR = "#0072B2"),
  Transplant = c(NS = "#D55E00", BR = "#009E73", FR = "#0072B2"))

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS6a_source_sym_NS_BR_DEGS_annotated.pdf",height=3,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames =F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
#pdf("heatmap_NSBR_source_sym.pdf",height=4,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()
####################################################################################################
#####trans NS versus FR
head(rldpvals)
p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.tNSFR<=p.val & !is.na(rld_data$padj.tNSFR),]
length(conds[,1])
#0

#source NS vs FR
colnames(rldpvals)
rld_data=as.data.frame(rldpvals[,c(1, 2, 3,  5, 6,7,8, 10, 12, 13,14, 15, 18, 21, 22,24,28,29,30,31,33,36,37,38,40, 41, 42)])
head(rld_data)
colnames(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.sNSFR<=p.val & !is.na(rld_data$padj.sNSFR),]
length(conds[,1])
#29
exp=conds[,1:25]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:26]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:26]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

colnames(exp)
my_sample_col <- data.frame(Source = c("FR", "FR","NS", "FR","FR", "NS", "FR", "NS", "NS", "NS", "NS", "FR","FR",   "FR", "NS","FR","FR","FR","FR", "NS", "NS", "NS", "FR", "NS", "NS"), Transplant = c("FR","BR","FR", "NS", "BR", "FR","BR", "NS","FR","FR", "NS", "FR","FR","BR","BR","FR","BR","BR","FR","BR","BR", "BR", "FR" ,"NS", "NS"))
row.names(my_sample_col)= colnames(anno)

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS6b_source_sym_NS_FR_DEGS_annotated.pdf",height=6,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames = F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
#pdf("heatmap_NSFR_source_sym.pdf",height=4,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()

####################################################################################################
#####trans BR versus FR
colnames(rldpvals)
rld_data=as.data.frame(rldpvals[,c(1, 2, 3,4, 6,7,8,11, 12, 13, 15,16, 18,19,20, 21, 22,23,24,25,26,28,29,30,31,32,33,34,35,36,37,39,51,52)])
head(rld_data)
colnames(rld_data)

p.val=0.10 # FDR cutoff
head(rld_data)
conds=rld_data[rld_data$padj.tBRFR<=p.val & !is.na(rld_data$padj.tBRFR),]
length(conds[,1])
#3
head(conds)
exp=conds[,1:32]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:33]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:33]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))
#no genes were annotated
# big heat map of all DEGs, no annotation too
#pdf("heatmap_NSFR_source_sym.pdf",height=4,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()

#source BR vs FR
head(rldpvals)
colnames(rldpvals)
rld_data=rldpvals[,c(1, 2, 4,5, 6,  8,9, 11, 15,16, 17,18, 19, 20, 21, 23,24, 25, 26, 27,28,29, 30,  32, 34, 35,37,39, 45, 46)]
head(rld_data)
colnames(rld_data)

p.val=0.10 # FDR cutoff
conds=rld_data[rld_data$padj.sBRFR<=p.val & !is.na(rld_data$padj.sBRFR),]
length(conds[,1])
#73
colnames(conds)
exp=conds[,1:28]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them

df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:29]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:29]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

colnames(exp)
my_sample_col <- data.frame(Source = c("FR","FR", "BR", "FR","FR", "FR", "BR","BR", "FR", "BR", "BR", "FR", "BR", "BR",  "FR", "BR", "FR", "BR","BR" , "BR","FR","FR", "FR","BR",  "BR","BR", "FR", "BR"), Transplant = c("FR","BR","FR","NS", "BR", "BR","NS","BR",  "FR","BR", "NS", "FR","FR","BR","BR","BR","FR","FR", "FR", "NS", "BR", "BR", "FR","FR","BR", "FR","FR","BR"))
row.names(my_sample_col)= colnames(anno)

# big heat map of all annotated genes
pdf("Figures/Supplemental_Figures/For_edits/FigS6c_source_BR_FR_DEGS_annotated_sym.pdf",height=8,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames = F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
#pdf("heatmap_BRFR_source_sym.pdf",height=5,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
#dev.off()

###########################################################################################
# VENN Diagram to include both up and down regulated genes in common for source
# install.packages('VennDiagram', dependencies=TRUE, repos='http://cran.us.r-project.org')
library(VennDiagram)
head(resSource1)
head(resSource2)
head(resSource3)

pSource1up=row.names(resSource1[resSource1$padj<0.1 & !is.na(resSource1$padj) & resSource1$log2FoldChange>0,])
length(pSource1up) #16
pSource1down=row.names(resSource2[resSource1$padj<0.1 & !is.na(resSource1$padj) & resSource1$log2FoldChange<0,])
length(pSource1down) #13
pSource2up=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj) & resSource2$log2FoldChange>0,])
length(pSource2up) #5
pSource2down=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj) & resSource2$log2FoldChange<0,])
length(pSource2down) #4
pSource3up=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj) & resSource3$log2FoldChange>0,])
length(pSource3up) #28
pSource3down=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj) & resSource3$log2FoldChange<0,])
length(pSource3down) #45

pSource1=row.names(resSource1[resSource1$padj<0.1 & !is.na(resSource1$padj),])
pSource2=row.names(resSource2[resSource2$padj<0.1 & !is.na(resSource2$padj),])
pSource3=row.names(resSource3[resSource3$padj<0.1 & !is.na(resSource3$padj),])

#UP
pdegs1_up=union(pSource1up,pSource2up)
pdegs1_up=union(pdegs1_up,pSource3up)
length(pdegs1_up)
#42

#DOWN
pdegs1_down=union(pSource1down,pSource2down)
pdegs1_down=union(pdegs1_down,pSource3down)
length(pdegs1_down)
#56

#ALL
pdegs1=union(pSource1,pSource2)
pdegs1=union(pdegs1,pSource3)
length(pdegs1)
#94

###do UP, DOWN, ALL
candidates=list("NSvsFRup"=pSource1up, "NSvsBRup"=pSource2up, "BRvsFRup"=pSource3up)
#quartz()
#prettyvenn=venn.diagram(
#  x = candidates,
#  filename=NULL,
#  col = "transparent",
#  fill = c("coral2", "forestgreen", "royalblue1"),
#  alpha = 0.5,
#  label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
#  cex = 2.5,
#  fontfamily = "sans",
#  fontface = "bold",
#  cat.default.pos = "text",
#  cat.col = c("darkred", "darkgreen", "blue4"),
#  cat.cex = 2.5,
#  cat.fontfamily = "sans",
#  cat.dist = c(0.08, 0.08, 0.03),
#  cat.pos = 1
#)
#grid.draw(prettyvenn)

# VENN Diagram to include both up and down regulated genes in common for transplant
head(resTrans1)
head(resTrans2)
head(resTrans3)

pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #0
pTrans1down=row.names(resTrans2[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #0
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #0
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #0
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #0
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #3

pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])

#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
length(pdegs1_up)
#0

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
length(pdegs1_down)
#3

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
length(pdegs1)
#3

###do UP, DOWN, ALL
candidates=list("NSvsFRup"=pTrans1up, "NSvsBRup"=pTrans2up, "BRvsFRup"=pTrans3up)
#prettyvenn=venn.diagram(
#  x = candidates,
#  filename=NULL,
#  col = "transparent",
#  fill = c("coral2", "forestgreen", "royalblue1"),
#  alpha = 0.5,
#  label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
#  cex = 2.5,
#  fontfamily = "sans",
#  fontface = "bold",
#  cat.default.pos = "text",
#  cat.col = c("darkred", "darkgreen", "blue4"),
#  cat.cex = 2.5,
#  cat.fontfamily = "sans",
#  cat.dist = c(0.08, 0.08, 0.03),
#  cat.pos = 1
#);
#quartz()
#grid.draw(prettyvenn)

# CANONICAL CORRESPONDENCE ANALYSIS OF GENE EXPRESSION based on raw counts
##packages used in plasticity analysis. need to restart R to avoid DLL issue
library("ggpubr")
library("tidyverse")
library("dplyr")
library("performance")
library("finalfit")
library("see")
library("vegan")
library("adegenet") 
library("WGCNA") # for labels2colors

cc <- read.table("Data/GeneExpression_data/sym_rtssid_counts_BM3.txt")
head(cc)
names(cc)=sub(".fastq.trim.sam.counts","",names(cc))
names(cc)

#remove outliers
cc$Ssid401 <- NULL
cc$Ssid447 <- NULL
cc$Ssid467 <- NULL
cc$Ssid498 <- NULL
ncol(cc)

mns=apply(cc,1,mean)
table(mns>1)
cc=cc[mns>1,]
cct=t(cc) # in vegan, species (genes) must be columns, and sites (samples) must be rows
head(cct)
names(cc)
length(cc)

rt=read.csv("Data/GeneExpression_data/RT_metadata_RTonly.csv")
rt2=rt[!(rt$sample %in% c("Ssid401", "Ssid447", "Ssid467", "Ssid498")), ]
conditions=data.frame(rt2$source, rt2$transplant)
nrow(conditions) 
names(conditions)=c("source", "trans")

library(vegan)
# normalizing to total counts
cct=decostand(cct,method="total")
# log-transform leaving zeroes as zeroes
cct=decostand(cct,method="log")

# canonical correspondence analysis (chi-square distances) - prepare to wait a while
ccaa=cca(cct~., conditions)
head(ccaa)
summary(ccaa)
# Accumulated constrained eigenvalues
# Importance of components:
#   CCA1    CCA2     CCA3      CCA4
# Eigenvalue            0.001485 0.00127 0.001151 0.0009023
# Proportion Explained  0.308819 0.26409 0.239408 0.1876829
# Cumulative Proportion 0.308819 0.57291 0.812317 1.0000000

pt_colors <- ifelse(conditions$source == 'FR', '#0072B2',
                    ifelse(conditions$source == 'NS', '#D55E00', '#009E73'))
pt_transcol <- ifelse(conditions$trans == 'FR', '#0072B2',
                      ifelse(conditions$trans == 'NS', '#D55E00', '#009E73'))
shape_source <- ifelse(conditions$source == 'FR', 15,
                       ifelse(conditions$source == 'NS', 16, 17))
shape_trans <- ifelse(conditions$trans == 'FR', 15,
                      ifelse(conditions$trans == 'NS', 16, 17))

pdf("Figures/Supplemental_Figures/For_edits/Fig3b_CCA_sym_trans_final.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites", type="n", ylab="CCA2 (26.4%)", xlab="CCA1 (30.9%)", xlim=range(-2.3,3), ylim=range(-2.0,2.7))
points(ccaa,choices=c(1,2),col=pt_colors, pch=shape_trans)
ordispider(ccaa,groups=conditions$trans,col="grey50")
ordiellipse(ccaa, groups=conditions$trans, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 1.5, cex = 0.8, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("grey","grey","grey"), pch = c(16,17,15),
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.13, 0.005), title="Transplant")
dev.off()

pdf("Figures/Supplemental_Figures/For_edits/Fig3a_CCA_sym_source_final.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites",type="n",  ylab="CCA2 (26.4%)", xlab="CCA1 (30.9%)", xlim=range(-2.3,3), ylim=range(-2.0,2.7))
points(ccaa,choices=c(1,2),col=pt_colors,pch=shape_trans)
ordispider(ccaa,groups=conditions$source,col="grey50")
ordiellipse(ccaa, groups=conditions$source, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 1.5, cex = 0.8, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("grey","grey","grey"), pch = c(16,17,15),
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.13, 0.005), title="Transplant")
dev.off()

#----------------------------
# ANOVA
anova(ccaa,by = "term") 
# Permutation test for cca under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = cct ~ source + trans, data = conditions)
# Df ChiSquare      F Pr(>F)    
# source    2  0.002519 1.2720  0.001 ***
#   trans     2  0.002288 1.1553  0.004 ** 
#   Residual 35  0.034661                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## -------------------------------------------------------- Colleen's additions below for plasticity:
set.seed(16)
#### I decided it will be cleaner if you source these functions from another script so you will need to do that now! Just save it to your working directory or add the path to it in the source function:
source("Code/CCAplast_function.R")

screeplot(ccaa)
## create dataframe of CCAs with experimental conditions
ccaa_full <- data.frame(conditions,
                        trt = paste(conditions$source, conditions$trans, sep = "_"),
                        scores(ccaa, choices = c(1:4))$sites)
# create dataframes for each source location
ccaa_fr <- ccaa_full %>% filter(trans == "FR")
ccaa_br <- ccaa_full %>% filter(trans == "BR")
ccaa_ns <- ccaa_full %>% filter(trans == "NS")

#### Calculate plasticity for each source location with custom function:
## Right now, these are calculating the plasticity based on all CCAs
# If you want to change that, you can specify it by unhashing 'num_cca' and specifying how many CCAs to include

plast_symb_FR <-  CCAplast(cca = ccaa_fr[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_fr[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "FR_FR", # control level of the treatment. If blank, a control mean per control level is assumed
                      keep = "yes") # change to 'yes' if you want to retain the control level values

plast_symb_BR <-  CCAplast(cca = ccaa_br[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_br[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "BR_BR", # control level of the treatment. If blank, a control mean per control level is assumed   
                      keep = "yes") # change to 'yes' if you want to retain the control level values

plast_symb_NS <-  CCAplast(cca = ccaa_ns[,-c(1:3)], # the CCA dataframe containing the CCA eigenvalues
                      data = ccaa_ns[,c(1:3)], # the condition/treatment data corresponding to samples
                      # sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_cca =  "2", # the number of CCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
                      control_col = "trt", # what the 'treatment' column is called
                      control_lvl = "NS_NS", # control level of the treatment. If blank, a control mean per control level is assumed
                      keep = "yes") # change to 'yes' if you want to retain the control level values

## combine all three plasticity dataframes into one
ccaa_symb_plast <- rbind(plast_symb_FR, plast_symb_BR, plast_symb_NS)
#table(ccaa_symb_plast$trt)
ccaa_symb_plast <- subset(ccaa_symb_plast, trt != "FR_NS") %>%  droplevels() # dropping the FR to NS level since it only has a single frag


##############
## Selecting best-fit model of GLM via AIC
# need to first filter out the single FR -> NS sample
# if redoing without outliers removed, do not need to do this ----->
#ccaa_plast <- ccaa_plast %>% filter(trt != "FR_NS") %>% droplevels()
## run a full model
plast_symb_mod <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_symb_plast)
#check_model(plast_symb_mod)

## data look fairly normal so proceeding with GLM here ##
## fit different model combos
plast_symb_mod <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_symb_plast) # best-fit model
plast_symb_mod2 <- glm(dist ~ source + trans, family = Gamma(link = "log"), data = ccaa_symb_plast)
plast_symb_mod3 <- glm(dist ~ trans, family = Gamma(link = "log"), data = ccaa_symb_plast)
plast_symb_mod4 <- glm(dist ~ source, family = Gamma(link = "log"), data = ccaa_symb_plast)
plast_symb_mod5 <- glm(dist ~ trt, family = Gamma(link = "log"), data = ccaa_symb_plast)

## check for best-fit model
compare_performance(plast_symb_mod, plast_symb_mod2, plast_symb_mod3, plast_symb_mod4, plast_symb_mod5)
# Name            | Model |    AIC |  AIC_wt |     BIC |  BIC_wt | Nagelkerke's R2 |  RMSE | Sigma
# ------------------------------------------------------------------------------------------------
# plast_symb_mod  |   glm | 31.742 |   0.500 |  46.714 |   0.500 |           0.863 | 0.273 | 0.368
# plast_symb_mod2 |   glm | 94.949 | < 0.001 | 104.930 | < 0.001 |           0.038 | 0.765 | 0.821
# plast_symb_mod3 |   glm | 91.827 | < 0.001 |  98.481 | < 0.001 |           0.011 | 0.763 | 0.807
# plast_symb_mod4 |   glm | 91.317 | < 0.001 |  97.971 | < 0.001 |           0.027 | 0.756 | 0.802
# plast_symb_mod5 |   glm | 31.742 |   0.500 |  46.714 |   0.500 |           0.863 | 0.273 | 0.368

## Best-fit GLMM with Gamma log link
plast_symb_glm <- glm(dist ~ source * trans, family = Gamma(link = "log"), data = ccaa_symb_plast)
summary(plast_symb_glm) # summary output

# Call:
#   glm(formula = dist ~ source * trans, family = Gamma(link = "log"), 
#       data = ccaa_symb_plast)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.20193  -0.18062   0.00546   0.16344   0.66548  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -1.2588     0.1355  -9.290 1.80e-10 ***
#   sourceFR           1.7762     0.1916   9.270 1.90e-10 ***
#   sourceNS           1.5797     0.2142   7.374 2.66e-08 ***
#   transFR            1.5676     0.1916   8.181 3.06e-09 ***
#   transNS            2.0990     0.2347   8.944 4.29e-10 ***
#   sourceFR:transFR  -2.9321     0.2710 -10.820 4.72e-12 ***
#   sourceNS:transFR  -1.0877     0.3030  -3.590  0.00112 ** 
#   sourceFR:transNS       NA         NA      NA       NA    
# sourceNS:transNS  -3.0780     0.3319  -9.274 1.88e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.110152)
# 
# Null deviance: 23.6194  on 38  degrees of freedom
# Residual deviance:  4.2095  on 31  degrees of freedom
# AIC: 31.742
# 
# Number of Fisher Scoring iterations: 5
### Bootstrapping with custom function (bootFUN)
bootFUN <- function(model, newdata) {
  if(class(model)[1] == "negbin"){
    nr <- nrow(model[["model"]]) # count number of rows in the model (# of observations)
    data <- model[["model"]] # pull data from model
  } else {
    nr <- nrow(model$data) # count number of rows in the model (# of observations)
    data <- model$data # pull data from model
  }
  data2 <- data.frame(data[sample(1:nr, replace = TRUE), ]) # random sample of numbers from 1 - nr (# of rows) to select data from
  update <- update(model, data = data2) # rerun the model using the 'new' dataset from the random smapling above
  predict(update, newdata, type = "response", level = 0.95, allow.new.levels = TRUE) # predicts response variable with model and updated data 
}
bootnum <- 1500 # this is the number of bootstrap replicates. I set it to 1500 here but we should select a number between 999 and 9999

# construct new dataframe and then run bootstrapping
newdata_symb_boot <- data.frame(coral = rownames(ccaa_symb_plast), trt = ccaa_symb_plast$trt, source = ccaa_symb_plast$source, trans = ccaa_symb_plast$trans, dist = ccaa_symb_plast$dist)
boot_symb <- replicate(bootnum, bootFUN(model = plast_symb_glm, newdata = newdata_symb_boot)) # run the bootFUN the number (bootnum) of times specified and save as matrix

# Calculate the mean, 95% lowerCI, and 95% upperCI from the boot matrix and add it to dataframe (but as a new name)
plast_symb_boot <- cbind(ccaa_symb_plast, as.data.frame(t(apply(boot_symb, 1, function(x) c(quantile(x, c(0.025, 0.5, 0.975)))))))
colnames(plast_symb_boot)[5:7] <- c("lowerci", "estimate", "upperci") # rename mean/CI columns


##############

# plot of plasticity (symbionts)
plast_symb_boot$source <- relevel(plast_symb_boot$source, ref = "NS")
plast_symb_boot$trans <- relevel(plast_symb_boot$trans, ref = "NS")
plast_symb_boot$trt <- paste(plast_symb_boot$source, plast_symb_boot$trans, sep = "_")
plast_symb_boot <- subset(plast_symb_boot, trt != "FR_NS") %>%  droplevels() # dropping the FR to NS level since it only has a single frag

ggplot(data = plast_symb_boot, aes(x = trans, colour = source)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 0, 0))), line = "none") +  
  geom_point(aes(y = dist), position = position_jitterdodge(dodge.width = 0.3), alpha = 0.4, size = 1, shape = 1) +
  geom_linerange(aes(ymin = lowerci, ymax = upperci, colour = source), size = 0.7, position = position_dodge(width = 0.3)) +
  geom_point(aes(y = estimate), position = position_dodge(width = 0.3), size = 2) +
  scale_color_manual("Source", values = c("#D55E00", "#009E73", "#0072B2")) +
  labs(x = "Transplant Location", y = "Plasticity")
ggsave("Figures/Supplemental_Figures/For_edits/Fig3c_CCA_sym_plasticity_CI.pdf", height = 3, width = 3, useDingbats = FALSE)

ccaa_symb_plast_sum <- ccaa_symb_plast %>%
  group_by(trans, source) %>%
  summarise(mean = mean(dist, na.rm = TRUE),
            sd = sd(dist))
#relevel
ccaa_symb_plast$source <- relevel(ccaa_symb_plast$source, ref = "NS")
ccaa_symb_plast$trans <- relevel(ccaa_symb_plast$trans, ref = "NS")
ccaa_symb_plast_sum$source <- relevel(ccaa_symb_plast_sum$source, ref = "NS")
ccaa_symb_plast_sum$trans <- relevel(ccaa_symb_plast_sum$trans, ref = "NS")

# plot of plasticity
#pdf("CCA_sym_plasticity_SD.pdf",height=3,width=3, useDingbats=FALSE)
ggplot(data = ccaa_symb_plast, aes(x = trans, colour = source)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 0, 0))), line = "none") +  
  geom_point(aes(y = dist), position =position_jitterdodge(dodge.width = 0.3), alpha = 0.4, size = 1, shape = 1) +
  geom_errorbar(data = ccaa_symb_plast_sum, aes(x = trans, ymin = mean - sd, ymax = mean + sd, colour = source), size = 0.7, width = 0.0, position = position_dodge(width = 0.2)) +
  geom_point(data = ccaa_symb_plast_sum, aes(y = mean), position = position_dodge(width = 0.2), size = 2) +
  scale_color_manual("Source", values = c("#D55E00", "#009E73", "#0072B2")) +
  labs(x = "Transplant Location", y = "Plasticity")
#dev.off()



####-------------------------------------------------------------------------------------------####
#### Colleen is adding in a heatmap of interesting symbiont genes from GO categories
# This is following what we did the in the ACER microplastics x OAW project: https://github.com/seabove7/Acer_OAW-Microplastics

##--------------get pvals 
# Source 1 = contrast=c("source","NS","FR")
valSource1 <- cbind(resSource1$pvalue, resSource1$padj)
head(valSource1)
colnames(valSource1) <- c("pval.Source1", "padj.Source1")
length(valSource1[,1])
table(complete.cases(valSource1))

# Source 2 = contrast=c("source","NS","BR")
valSource2 <- cbind(resSource2$pvalue, resSource2$padj)
head(valSource2)
colnames(valSource2) <- c("pval.Source2", "padj.Source2")
length(valSource2[,1])
table(complete.cases(valSource2))

# Source 3 = contrast=c("source","BR","FR")
valSource3 <- cbind(resSource3$pvalue, resSource3$padj)
head(valSource3)
colnames(valSource3) <- c("pval.Source3", "padj.Source3")
length(valSource3[,1])
table(complete.cases(valSource3))


######-------------make rlogdata and pvals table
rlog_symb <- rlogTransformation(dds_symb, blind = TRUE) 
rld_symb <- assay(rlog_symb)
head(rld_symb)
colnames(rld_symb) <- paste(conditions$source, conditions$trans, sep = "_")
head(rld_symb)
length(rld_symb[,1])

rldpvals <- cbind(rld_symb, valSource1, valSource2, valSource3)
head(rldpvals)
dim(rldpvals)
# 4750   46
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 3315  1435 

write.csv(rldpvals, "Data/GeneExpression_data/Belize_RT_RLDandPVALS.csv", quote=F)




library(dplyr)
library(stringr)
iso2go <- read.table("GO_enrichment_Symb/sym_feb_iso2go.tab", fill = TRUE) %>%
  dplyr::rename("GO_ID" = "V2") 
  
head(iso2go)

sym_source_BR_FR_res <- read.table("Data/GeneExpression_data/sym_RTsource_BRvsFR.txt")
head(sym_source_BR_FR_res)


rldpval <- readr::read_csv("Data/GeneExpression_data/Belize_RT_RLDandPVALS.csv") %>%
  select(gene = 1, everything())
head(rldpval)

col0 <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

##pulling all of the immunity terms and making a data frame
symb_source <- iso2go %>%
  mutate(gene = iso2go$V1) %>%
  left_join(rldpval)
head(symb_source)

#immuno=symb_source[,c(2,3,6,8,11,12,13,15,17,19,20,26,27)]
#head(immuno)
row.names(symb_source) <- symb_source$gene
symb_source <- symb_source[,-1]
head(symb_source)

##gene table
gg <- read.table("GO_enrichment_Symb/sym_feb_iso2gene.tab",sep="\t")
head(gg)

p.val <- 0.10 # raw pvalue for GO enriched
conds <- symb_source[symb_source$pval.Source1<=p.val & !is.na(symb_source$pval.Source1),]
length(conds[,1])
#387
head(conds)

############## COLLEEN STOPPED HERE ##############


exp <- conds[,1:10]
head(exp)
means <- apply(exp,1,mean) # means of rows
explc <- exp - means # subtracting them
head(explc)
head(gg)

library(tidyverse)

df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:11]

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:11]
head(anno)

#dataframe of the samples
colnames(anno)

##write out and edit manually a few gene names that have weird descriptions
write.csv(anno, "immunity_info_raw.csv", quote=TRUE)
anno2 <- read.csv("immunity_info_EDIT.csv", row.names=1) # Sarah performed some modifications to the excel file before uploading it here

#cbPalette2 <- c("darkorange","firebrick2", "firebrick4")
my_sample_col <- data.frame(treatment = c("OAW+MP","AMB","AMB","OAW+MP","AMB","OAW+MP","OAW+MP", "OAW+MP","AMB","AMB"))
row.names(my_sample_col)= colnames(anno2)
my_colour = list(treatment = c(`OAW+MP` = cbPalette[4], AMB = cbPalette[1]))

# big heat map of all annotated genes
library(pheatmap)
pdf("Figures/Figure4_OAW_AMB_immunity.pdf", height = 8.2, width = 9, onefile = F)
pheatmap(anno2, cluster_cols = TRUE, scale = "row", color = col0, annotation_col = my_sample_col, annotation_colors = my_colour, show_rownames = TRUE, show_colnames = FALSE, border_color = "NA")
dev.off()

setwd("~/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
bams=read.table("ssid_bams.txt")[,1] # list of bam files
goods=c(1:length(bams))
head(goods)

#--------------------
# loading individual to population correspondences
i2p=read.table("ssid_bams_indiv.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
head(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

ma = as.matrix(read.table("ssid.ibsMat"))

ma=ma[goods,goods]
head(ma)
length(ma)
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf('Host_Dendro_AllData.pdf', width = 8, height = 5)
plot(hc,cex=0.7) 
dev.off()
###Looks like there are some outliers and one pair might be clonal

# performing PCoA and CAP
library(vegan)
head(site)
pp0=capscale(ma~1)
pp=capscale(ma~site)

# significance of by-site divergence
adonis2(ma~site)
# Terms added sequentially (first to last)
# 
# #           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# site       2   0.07132 0.03566 0.94641 0.04413  0.611
# Residuals 41   1.54487 0.03768         0.95587       
# Total     43   1.61619                 1.00000  

# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
library(adegenet) # for transp()
cmd=pp
conds=as.data.frame(site)
conds$col <- ifelse(conds$site == "NS", "#D55E00", ifelse(conds$site == "BR", "#009E73", "#0072B2"))
pdf('capscale_Host_PCA_AllData.pdf', width = 6, height = 6)
plot(cmd,choices=axes2plot,pch=19,col=conds$col,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=conds$col)
ordispider(cmd,choices=axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices=axes2plot,groups= conds$site, draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)
dev.off()

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=conds$col)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=8,cex=0.7)
##we see that 470, 471, 472, were outliers based on PCA
##removed samples and re-run angsd pipeline. 

##################################################################################################################################################################################################################
###PART II host with outliers removed
setwd("~daviessw/Library/CloudStorage/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
bams=read.table("bams_remove.txt")[,1] # list of bam files
goods=c(1:length(bams))
head(goods)

#--------------------
# loading individual to population correspondences
i2p=read.table("ssid_bams_remove_indiv.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
head(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
ma = as.matrix(read.table("ssid_remove.ibsMat"))

ma=ma[goods,goods]
head(ma)
length(ma)
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf('Host_Dendro_removed.pdf', width = 8, height = 5)
plot(hc,cex=0.7) 
dev.off()

# performing PCoA and CAP
library(vegan)
head(site)
pp0=capscale(ma~1)
pp=capscale(ma~site)

# significance of by-site divergence
adonis2(ma~site)
          # Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
# site       2   0.07526 0.037630   1.025 0.0525  0.296
# Residuals 37   1.35835 0.036712         0.9475       
# Total     39   1.43361                  1.0000     

axes2plot=c(1,2)  
library(adegenet) # for transp()
cmd=pp
conds=as.data.frame(site)
conds$col <- ifelse(conds$site == "NS", "#D55E00", ifelse(conds$site == "BR", "#009E73", "#0072B2"))
plot(cmd,choices=axes2plot,pch=19,col=conds$col,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=conds$col)
ordispider(cmd,choices=axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices=axes2plot,groups= conds$site, draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=conds$col)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=8,cex=0.7)
##we clearly see that 413 and 476 are being identified as high liklihood of sharing common ancestry
##removed 413 and re-run angsd pipeline with 39 samples 

###PART III host with potential clone 413) removed
setwd("~daviessw/Library/CloudStorage/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
bams=read.table("bams_remove2.txt")[,1] # list of bam files
goods=c(1:length(bams))
head(goods)

#--------------------
# loading individual to population correspondences
i2p=read.table("ssid_bams_remove2_indiv.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
head(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
ma = as.matrix(read.table("ssid_remove2.ibsMat"))

ma=ma[goods,goods]
head(ma)
length(ma)
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf('Host_Dendro_removed2.pdf', width = 8, height = 5)
plot(hc,cex=0.7) 
dev.off()

# performing PCoA and CAP
library(vegan)
head(site)
pp0=capscale(ma~1)
pp=capscale(ma~site)

pp=cca(ma~site)

# significance of by-site divergence
adonis2(ma~site)
          # Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# site       2   0.07037 0.035187 0.96608 0.05094  0.792
# Residuals 36   1.31119 0.036422         0.94906       
# Total     38   1.38156                  1.00000      

axes2plot=c(1,2)  
library(adegenet) # for transp()
cmd=pp
conds=as.data.frame(site)
conds$col <- ifelse(conds$site == "NS", "#D55E00", ifelse(conds$site == "BR", "#009E73", "#0072B2"))
plot(cmd,choices=axes2plot,pch=19,col=conds$col,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=conds$col)
ordispider(cmd,choices=axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices=axes2plot,groups= conds$site, draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)

# unscaled, to identify outliers
pdf('Host_PCA_removed2.pdf', width = 5, height = 5)
plot(cmd$CA$u[,axes2plot],pch=19,col=conds$col)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)
dev.off()

###try CCA like we did for Gene Expression for response to reviewers
# canonical correspondence analysis 
bams=read.table("bams_remove2.txt")[,1] # list of bam files
goods=c(1:length(bams))
head(goods)

i2p=read.table("ssid_bams_remove2_indiv.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
head(i2p)
nrow(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
co = as.matrix(read.table("ssid_remove2.ibsMat"))
head(co)
nrow(co)
head(conds)
site=data.frame(site)
ccaa=cca(co~., site)
head(ccaa)
summary(ccaa)

#Partitioning of scaled Chi-square:
#Inertia Proportion
#Total         0.026764    1.00000
#Constrained   0.001359    0.05079
#Unconstrained 0.025405    0.94921

#Importance of components:
#                       CCA1      CCA2      
#Eigenvalue            0.00071 0.0006493
#Proportion Explained  0.02653 0.0242599 
#Cumulative Proportion 0.02653 0.0507897

#CCA1 var= 0.02653
#CCA2 var=0.02426

anova(ccaa,by = "term") 
#Model: cca(formula = co ~ conds.site, data = site)
#Df ChiSquare      F Pr(>F)
#conds.site  2 0.0013593 0.9631  0.809
#Residual   36 0.0254046

conds=as.data.frame(site)
conds$col <- ifelse(conds$site == "NS", "#D55E00", ifelse(conds$site == "BR", "#009E73", "#0072B2"))

pdf("CCA_Host_SNP.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites", type="n", ylab="CCA2 (2.4%)", xlab="CCA1 (2.7%)", xlim=range(-2,2.2), ylim=range(-2,2))
points(ccaa,choices=c(1,2),col=conds$col, pch = 16)
ordispider(ccaa,groups=conds$site,col="grey50")
ordiellipse(ccaa, groups=site$site, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 2, cex = 1, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
dev.off()


###Now the algal symbiont###
##PART I symbiont all 44 samples
setwd("~daviessw/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
bams=read.table("ssid_bams.txt")[,1] # list of bam files
goods=c(1:length(bams))
head(goods)

#--------------------
# loading individual to population correspondences
i2p=read.table("ssid_bams_indiv.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
head(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
ma = as.matrix(read.table("sym.ibsMat"))

ma=ma[goods,goods]
head(ma)
length(ma)
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf('Sym_Dendro.pdf', width = 8, height = 5)
plot(hc,cex=0.7) 
dev.off()

# performing PCoA and CAP
library(vegan)
head(site)
pp0=capscale(ma~1)
pp=capscale(ma~site)

# significance of by-site divergence
adonis2(ma~site)
          # Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# site       2   0.11655 0.058275  1.1688 0.05394  0.096 .
# Residuals 41   2.04420 0.049858         0.94606         
# Total     43   2.16075                  1.00000        

axes2plot=c(1,2)  
library(adegenet) # for transp()
cmd=pp
conds=as.data.frame(site)
conds$col <- ifelse(conds$site == "NS", "#D55E00", ifelse(conds$site == "BR", "#009E73", "#0072B2"))
plot(cmd,choices=axes2plot,pch=19,col=conds$col,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=conds$col)
ordispider(cmd,choices=axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices=axes2plot,groups= conds$site, draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)

# unscaled, to identify outliers
pdf('sym_PCA.pdf', width = 5, height = 5)
plot(cmd$CA$u[,axes2plot],pch=19,col=conds$col)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=c("#009E73", "#0072B2","#D55E00"),label=T)
dev.off()

###try CCA like we did for Gene Expression for response to reviewers
# canonical correspondence analysis 
head(i2p)
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]
head(i2p)
nrow(i2p)
co = as.matrix(read.table("sym.ibsMat"))
head(co)
head(conds)
site=data.frame(conds$site)
site
ccaa=cca(co~., site)
head(ccaa)
summary(ccaa)
#Partitioning of scaled Chi-square:
#  Inertia Proportion
#Total         0.026904    1.00000
#Constrained   0.001478    0.05492
#Unconstrained 0.025427    0.94508

#Eigenvalues, and their contribution to the scaled Chi-square 

#Importance of components:
 #                         CCA1      CCA2      
#Eigenvalue            0.0007958 0.0006818 
#Proportion Explained  0.0295784 0.0253409 
#Cumulative Proportion 0.0295784 0.0549193

#CCA1 var= 0.0295784
#CCA2 var=0.0253409

anova(ccaa,by = "term") 
#Model: cca(formula = co ~ conds.site, data = site)
#Df ChiSquare      F Pr(>F)  
#conds.site  2 0.0014776 1.1913  0.079 .
#Residual   41 0.0254266 

pdf("CCA_Sym_SNP.pdf",height=5,width=5, useDingbats=FALSE)
plot(ccaa,choices=c(1,2),display="sites", type="n", ylab="CCA2 (2.5%)", xlab="CCA1 (3.0%)", xlim=range(-2.5,2.5), ylim=range(-3,3))
points(ccaa,choices=c(1,2),col=conds$col, pch = 16)
ordispider(ccaa,groups=conds$site,col="grey50")
ordiellipse(ccaa, groups=conds$site, kind="sd", lwd=2, col=c("#009E73", "#0072B2", "#D55E00"), alpha=0.5, label=F)
legend("bottomleft", legend = c("NS", "BR", "FR"), col = c("#D55E00","#009E73", "#0072B2"), pch = 16,  bty = "n", pt.cex = 2, cex = 1, text.col = "black", horiz = F , inset = c(0.005, 0.005), title="Source")
dev.off()


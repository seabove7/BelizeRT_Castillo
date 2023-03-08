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
adonis(ma~site)
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
setwd("~daviessw/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
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
adonis(ma~site)
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
setwd("~daviessw/Dropbox/UNC/RT_Paper/2022_Final_RT/SNP/2022_SNP")
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

# significance of by-site divergence
adonis(ma~site)
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
adonis(ma~site)
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

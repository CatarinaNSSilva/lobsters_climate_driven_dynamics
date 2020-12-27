library(adegenet)
library(pegas) # funny thing reading files with numbers in beginng/end of file name

# vcf  JTR_JPA batch 1 ####
vcf_JTR_JPA = read.vcf("./output/batch_1_JTR_JPA.vcf", which.loci = 1:1e6)
dim(vcf_JTR_JPA) #  98  indv, 18879 SNPs

genind_JTR_JPA <- loci2genind(vcf_JTR_JPA)
summary_JTR_JPA = summary(genind_JTR_JPA)

str(summary_JTR_JPA)
summary_JTR_JPA$n.by.pop
mean(summary_JTR_JPA$Hobs) # 0.1739588
mean(summary_JTR_JPA$Hexp) # 0.2250516
min(summary_JTR_JPA$loc.n.all)
max(summary_JTR_JPA$loc.n.all)
mean(summary_JTR_JPA$NA.perc) # 28.86584

df_JTR_JPA = genind2df(genind_JTR_JPA, sep = "/")
df_JTR_JPA[1:15,1:5]

write.csv(df_JTR_JPA, file = "./output/df_JTR_JPA.csv", row.names = T)
df_JTR_JPA_SNPs = read.csv("./output/df_JTR_JPA.csv", header = F, row.names = 1)
df_JTR_JPA_SNPs[1:15,1:5]
df_JTR_JPA_SNPs = df_JTR_JPA_SNPs[1,]
write.csv(df_JTR_JPA_SNPs, file = "./output/df_JTR_JPA_SNPs.csv", row.names = F)

# strata 2 pops ####
JTR_JPA_strata <- read.table("./input/strata_JTR_JPA.txt", header=T)
head(JTR_JPA_strata)
dim(JTR_JPA_strata)

strata(genind_JTR_JPA) <- JTR_JPA_strata
genind_JTR_JPA

library("mmod")
setPop(genind_JTR_JPA) = ~pop
str(genind_JTR_JPA)


summary_JTR_JPA = summary(genind_JTR_JPA) # 98 indv, 2 pops
str(summary_JTR_JPA)


# strata 5 pops #### 

JTR_JPA_strata_5pops <- read.table("./input/strata_JTR_JPA_5pops.txt", header=T)
head(JTR_JPA_strata_5pops)
dim(JTR_JPA_strata_5pops)

strata(genind_JTR_JPA) <- JTR_JPA_strata_5pops
genind_JTR_JPA

library("mmod")
setPop(genind_JTR_JPA) = ~pop
str(genind_JTR_JPA)

summary_JTR_JPA_5pops = summary(genind_JTR_JPA) # 98 indv, 5 pops
summary_JTR_JPA_5pops$n.by.pop




# HO HE ####
# adegenet #
genind_JTR_JPA$pop
# Levels: JPA_SP JPA_AM JTR_TC-II JTR_GI JTR_NI

JPA_SP_genind = seppop(genind_JTR_JPA)$JPA_SP
summary_JPA_SP_genind = summary(JPA_SP_genind)
str(summary_JPA_SP_genind)
mean(summary_JPA_SP_genind$Hobs, na.rm=T)
mean(summary_JPA_SP_genind$Hexp, na.rm=T)

JPA_AM_genind = seppop(genind_JTR_JPA)$JPA_AM
summary_JPA_AM_genind = summary(JPA_AM_genind)
str(summary_JPA_AM_genind)
mean(summary_JPA_AM_genind$Hobs, na.rm=T)
mean(summary_JPA_AM_genind$Hexp, na.rm=T)

JTR_TC_genind = seppop(genind_JTR_JPA)$JTR_TC
summary_JTR_TC_genind = summary(JTR_TC_genind)
str(summary_JTR_TC_genind)
mean(summary_JTR_TC_genind$Hobs, na.rm=T)
mean(summary_JTR_TC_genind$Hexp, na.rm=T)

JTR_GI_genind = seppop(genind_JTR_JPA)$JTR_GI
summary_JTR_GI_genind = summary(JTR_GI_genind)
str(summary_JTR_GI_genind)
mean(summary_JTR_GI_genind$Hobs, na.rm=T)
mean(summary_JTR_GI_genind$Hexp, na.rm=T)

JTR_NI_genind = seppop(genind_JTR_JPA)$JTR_NI
summary_JTR_NI_genind = summary(JTR_NI_genind)
str(summary_JTR_NI_genind)
mean(summary_JTR_NI_genind$Hobs, na.rm=T)
mean(summary_JTR_NI_genind$Hexp, na.rm=T)



n.pop <- seppop(genind_JTR_JPA)
mySamp <- lapply(n.pop, function(x) mean(summary(x)$Hobs))

mean(summary(n.pop$JPA_SP))
mean(summary_JTR_JPA$Hobs)


mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 
# mean.hobs[is.nan(mean.hobs)] <- NA 
barplot(mean.hobs) 

mean.Hexp <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hexp))) 
# mean.Hexp[is.nan(mean.Hexp)] <- NA 






# DAPC ####

grp_JTR_JPA <- find.clusters(genind_JTR_JPA, max.n.clust=20) # retain 100 PCs, 2 clusters
dapc_JTR_JPA <- dapc(genind_JTR_JPA, grp_JTR_JPA$grp_JTR_JPA)
# 98 PCs,  $var = 1, a-score = 0
# retain 1 PCs (optimal),  $var = 0.0253541,  a-score = 0.4954248
# 2 PCs, $var = 0.04010284,  a-score = 0.4901552
# 4 PCs, $var = 0.07375324,  a-score = 0.3556437
# 5 PCs, $var = 0.08892232,  a-score = 0.3384488

head(dapc_JTR_JPA) # 


# a-score  ##

temp <- a.score(dapc_JTR_JPA)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score
temp$mean # 


png("./figs/a-score_optimisations_JTR_JPA.png", width=50, height=30, units="cm",res=155,pointsize = 36)
temp <- optim.a.score(dapc_JTR_JPA) # 1 PC optimal
dev.off()

myCol5 <- c("blue", "blue", "green", "green", "green") 

png("./figs/PCA_JTR_JPA_1PCs.png", width=30, height=20, units="cm",res=155,pointsize = 24)
scatter(dapc_JTR_JPA, posi.da="bottomright", bg="white", pch=17:22)
dev.off()

png("./figs/DAPC_JTR_JPA_4PCs.png", width=30, height=20, units="cm",res=155,pointsize = 18)
scatter(dapc_JTR_JPA, scree.da=FALSE, bg="white", pch=20, cell=1.5, cstar=0, 
        cex=3,clab=0, leg=T, txt.leg=paste(c("JPA_SP", "JPA_AM", "JTR_TC-II", "JTR_GI", "JTR_NI")), 
        posi.leg = "bottomleft", col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()


#  memb prob ####


png("./figs/memb_prob_JTR_JPA_1Pcs.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()



#  Ar ####

library(hierfstat)

JTR_JPA_hierfstat = genind2hierfstat(genind_JTR_JPA)
dim(JTR_JPA_hierfstat)
JTR_JPA_hierfstat[1:5,1:5]

allelic.richness_JTR_JPA = allelic.richness(JTR_JPA_hierfstat)
str(allelic.richness_JTR_JPA)
dim(allelic.richness_JTR_JPA$Ar)

AR_matrix_JTR_JPA = as.data.frame(allelic.richness_JTR_JPA$Ar)
str(AR_matrix_JTR_JPA)
dim(AR_matrix_JTR_JPA)
write.csv(AR_matrix_JTR_JPA, file = "./output/AR_matrix_JTR_JPA.csv")
# JPA = 1.856755		JTR = 1.847517

# to update #
JTR_JPA_hierfstat_5pops = genind2hierfstat(genind_JTR_JPA)
dim(JTR_JPA_hierfstat_5pops)
JTR_JPA_hierfstat_5pops[1:5,1:5]

allelic.richness_JTR_JPA = allelic.richness(JTR_JPA_hierfstat_5pops)
str(allelic.richness_JTR_JPA)
dim(allelic.richness_JTR_JPA$Ar)

AR_matrix_JTR_JPA = as.data.frame(allelic.richness_JTR_JPA$Ar)
str(AR_matrix_JTR_JPA)
dim(AR_matrix_JTR_JPA)
write.csv(AR_matrix_JTR_JPA, file = "./output/AR_matrix_JTR_JPA_5pops.csv")




# PCA ####

write.csv(JTR_JPA_hierfstat, file = "./output/JTR_JPA_hierfstat.csv", row.names = F)
# replace "NAs" for "0"
JTR_JPA_hierfstat_PCA = read.csv("./output/JTR_JPA_hierfstat_PCA_spp.csv", header = T)

JTR_JPA_hierfstat_PCA = JTR_JPA_hierfstat_PCA[,2:18880]
JTR_JPA_hierfstat_PCA[1:5,1:5]
dim(JTR_JPA_hierfstat_PCA)

JTR_JPA_indv_list = read.csv("./input/JTR_JPA_indv_list.txt", header = F)
JTR_JPA_indv_list = as.matrix(JTR_JPA_indv_list)
head(JTR_JPA_indv_list) 
JTR_JPA_hierfstat_PCA$id = JTR_JPA_indv_list
dim(JTR_JPA_hierfstat_PCA)
JTR_JPA_hierfstat_PCA[1:5,18875:18880] 
JTR_JPA_hierfstat_PCA = JTR_JPA_hierfstat_PCA[,c(18880,1:18879)]
str(JTR_JPA_hierfstat_PCA)
JTR_JPA_hierfstat_PCA[1:5,1:5]

write.csv(JTR_JPA_hierfstat_PCA, file = "./output/JTR_JPA_hierfstat_PCA.csv", row.names = F)
JTR_JPA_hierfstat_PCA = read.csv("./output/JTR_JPA_hierfstat_PCA.csv", header = T, row.names = 1)
JTR_JPA_hierfstat_PCA[1:5,1:5]



# row.names(JTR_JPA_hierfstat_PCA) = JTR_JPA_hierfstat_PCA[,1]

pca = prcomp(JTR_JPA_hierfstat_PCA)
summary(pca)
# PC1
# Standard deviation     903.3350
# Proportion of Variance   0.2874
# PC2
# Standard deviation    388.77997
# Proportion of Variance  0.05324



png("./figs/pca_JTR_JPA.png", width=40, height=30, units="cm",res=155,pointsize = 26)
plot(pca$x, pch=20, col="blue") #  xlim=c(-11,9), ylim=c(-4.5,4.5)
text(pca$x, rownames(pca$x), pos=3, cex=0.6)
dev.off()

strata = read.csv("./input/strata_JTR_JPA.txt")
library(ggfortify)
png("./figs/pca_JTR_JPA.png", width=15, height=8, units="cm",res=155,pointsize = 20)
autoplot(pca, data = strata, colour = 'pop'  ) + theme_bw() + scale_color_manual(values=c("blue", "green"))
dev.off()


strata_5pops = read.csv("./input/strata_JTR_JPA_5pops.txt")
library(ggfortify)
png("./figs/pca_JTR_JPA_5pops.png", width=15, height=8, units="cm",res=155,pointsize = 20)
autoplot(pca, data = strata_5pops, colour = 'pop'  ) + theme_bw() + scale_color_manual(values=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()



# vcf  outliers ####

system2("vcftools", 
        args = c("--vcf", "./output/batch_1_JTR_JPA.vcf",
                 "--stdout", "--recode",
                 "--recode-INFO-all", "--snps", "./output/outliers_Bayescan.txt",
                 ">", "./output/outliers_Bayescan.vcf"))


system2("vcftools", 
        args = c("--vcf", "./output/batch_1_JTR_JPA.vcf",
                 "--stdout", "--recode",
                 "--recode-INFO-all", "--snps", "./output/outliers_PCAdapt_Bayescan.txt",
                 ">", "./output/outliers_PCAdapt_Bayescan.vcf"))


vcf_JTR_JPA_outliers = read.vcf("./output/outliers_PCAdapt_Bayescan.vcf", which.loci = 1:1e6)
dim(vcf_JTR_JPA_outliers) #  98  indv, 1623 SNPs

genind_JTR_JPA_outliers <- loci2genind(vcf_JTR_JPA_outliers)
summary_JTR_JPA_outliers = summary(genind_JTR_JPA_outliers)

str(summary_JTR_JPA_outliers)
summary_JTR_JPA_outliers$n.by.pop
mean(summary_JTR_JPA_outliers$Hobs) # 0.1668602
mean(summary_JTR_JPA_outliers$Hexp) # 0.2124203
min(summary_JTR_JPA_outliers$loc.n.all)
max(summary_JTR_JPA_outliers$loc.n.all)
mean(summary_JTR_JPA_outliers$NA.perc) # 15.52052


# strata 5 pops #### 

JTR_JPA_strata_5pops <- read.table("./input/strata_JTR_JPA_5pops.txt", header=T)
head(JTR_JPA_strata_5pops)
dim(JTR_JPA_strata_5pops)

strata(genind_JTR_JPA_outliers) <- JTR_JPA_strata_5pops
genind_JTR_JPA_outliers

library("mmod")
setPop(genind_JTR_JPA_outliers) = ~pop
str(genind_JTR_JPA_outliers)

summary_JTR_JPA_5pops = summary(genind_JTR_JPA_outliers) # 98 indv, 5 pops
summary_JTR_JPA_5pops$n.by.pop


# DAPC ####

grp_JTR_JPA_outliers <- find.clusters(genind_JTR_JPA_outliers, max.n.clust=20) # retain 100 PCs, 4 clusters
dapc_JTR_JPA_outliers <- dapc(genind_JTR_JPA_outliers, grp_JTR_JPA_outliers$grp_JTR_JPA_outliers)
# 98/3= 33 PCs,       $var = 0.6077057, a-score = 
#  13 PCs (optimal),  $var = 0.3747567,  a-score = 0.3626128
# 2 PCs,              $var = 0.1702787,  a-score = 0.3299524
# 1 PC,               $var = 0.0975413,  a-score = 0.2888061
head(dapc_JTR_JPA_outliers) # 


# a-score  ##

temp <- a.score(dapc_JTR_JPA_outliers)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score
temp$mean # 


png("./figs/a-score_optimisations_JTR_JPA_outliers.png", width=50, height=30, units="cm",res=155,pointsize = 36)
temp <- optim.a.score(dapc_JTR_JPA_outliers) # 13 PC optimal
dev.off()

myCol5 <- c("darkblue", "lightblue", "green", "darkgreen", "lightgreen") 

png("./figs/PCA_JTR_JPA_1PCs.png", width=30, height=20, units="cm",res=155,pointsize = 24)
scatter(dapc_JTR_JPA_outliers, posi.da="bottomright", bg="white", pch=17:22)
dev.off()

png("./figs/DAPC_JTR_JPA_outliers_13PCs.png", width=30, height=20, units="cm",res=155,pointsize = 18)
scatter(dapc_JTR_JPA_outliers, scree.da=FALSE, bg="white", pch=20, cell=1.5, cstar=0, 
        cex=3,clab=0, leg=T, txt.leg=paste(c("JPA_SP", "JPA_AM", "JTR_TC-II", "JTR_GI", "JTR_NI")), 
        posi.leg = "bottomleft", col=myCol5)
dev.off()


#  memb prob ####


png("./figs/memb_prob_JTR_JPA_outliers_1Pcs_spp.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_outliers, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()



#  Ar outliers ####

library(hierfstat)

JTR_JPA_outliers_hierfstat = genind2hierfstat(genind_JTR_JPA_outliers)
dim(JTR_JPA_outliers_hierfstat)
JTR_JPA_outliers_hierfstat[1:5,1:5]

allelic.richness_JTR_JPA_outliers = allelic.richness(JTR_JPA_outliers_hierfstat)
str(allelic.richness_JTR_JPA_outliers)
dim(allelic.richness_JTR_JPA_outliers$Ar)

AR_matrix_JTR_JPA_outliers = as.data.frame(allelic.richness_JTR_JPA_outliers$Ar)
str(AR_matrix_JTR_JPA_outliers)
dim(AR_matrix_JTR_JPA_outliers)
write.csv(AR_matrix_JTR_JPA_outliers, file = "./output/AR_matrix_JTR_JPA_outliers.csv")
# 1.32	1.31	1.42	1.40	1.60

# to update #
JTR_JPA_hierfstat_5pops = genind2hierfstat(genind_JTR_JPA)
dim(JTR_JPA_hierfstat_5pops)
JTR_JPA_hierfstat_5pops[1:5,1:5]

allelic.richness_JTR_JPA = allelic.richness(JTR_JPA_hierfstat_5pops)
str(allelic.richness_JTR_JPA)
dim(allelic.richness_JTR_JPA$Ar)

AR_matrix_JTR_JPA = as.data.frame(allelic.richness_JTR_JPA$Ar)
str(AR_matrix_JTR_JPA)
dim(AR_matrix_JTR_JPA)
write.csv(AR_matrix_JTR_JPA, file = "./output/AR_matrix_JTR_JPA_5pops.csv")




# PCA outliers####

write.csv(JTR_JPA_outliers_hierfstat, file = "./output/JTR_JPA_outliers_hierfstat.csv", row.names = F)
# replace "NAs" for "0"
JTR_JPA_outliers_hierfstat_PCA = read.csv("./output/JTR_JPA_outliers_hierfstat_PCA_spp.csv", header = T)

JTR_JPA_outliers_hierfstat_PCA = JTR_JPA_outliers_hierfstat_PCA[,2:1624]
JTR_JPA_outliers_hierfstat_PCA[1:5,1:5]
dim(JTR_JPA_outliers_hierfstat_PCA) # 98 1623

JTR_JPA_indv_list = read.csv("./input/JTR_JPA_indv_list.txt", header = F)
JTR_JPA_indv_list = as.matrix(JTR_JPA_indv_list)
head(JTR_JPA_indv_list) 
JTR_JPA_outliers_hierfstat_PCA$id = JTR_JPA_indv_list
dim(JTR_JPA_outliers_hierfstat_PCA)
JTR_JPA_outliers_hierfstat_PCA[1:5,1620:1624] 
JTR_JPA_outliers_hierfstat_PCA = JTR_JPA_outliers_hierfstat_PCA[,c(1624,1:1623)]
str(JTR_JPA_outliers_hierfstat_PCA)
JTR_JPA_outliers_hierfstat_PCA[1:5,1:5]

write.csv(JTR_JPA_outliers_hierfstat_PCA, file = "./output/JTR_JPA_outliers_hierfstat_PCA.csv", row.names = F)
JTR_JPA_outliers_hierfstat_PCA = read.csv("./output/JTR_JPA_outliers_hierfstat_PCA.csv", header = T, row.names = 1)
JTR_JPA_outliers_hierfstat_PCA[1:5,1:5]



# row.names(JTR_JPA_outliers_hierfstat_PCA) = JTR_JPA_outliers_hierfstat_PCA[,1]

pca_outliers = prcomp(JTR_JPA_outliers_hierfstat_PCA)
summary(pca_outliers)
# PC1
# Standard deviation     149.4119
# Proportion of Variance   0.1319
# PC2
# Standard deviation    110.16123
# Proportion of Variance  0.07169



png("./figs/pca_outliers_JTR_JPA.png", width=40, height=30, units="cm",res=155,pointsize = 26)
plot(pca_outliers$x, pch=20, col="blue") #  xlim=c(-11,9), ylim=c(-4.5,4.5)
text(pca_outliers$x, rownames(pca_outliers$x), pos=3, cex=0.6)
dev.off()

strata = read.csv("./input/strata_JTR_JPA.txt")
library(ggfortify)
png("./figs/pca_outliers_JTR_JPA.png", width=15, height=8, units="cm",res=155,pointsize = 20)
autoplot(pca_outliers, data = strata, colour = 'pop'  ) + theme_bw() + scale_color_manual(values=c("blue", "green"))
dev.off()


strata_5pops = read.csv("./input/strata_JTR_JPA_5pops.txt")
library(ggfortify)
png("./figs/pca_outliers_JTR_JPA_5pops.png", width=15, height=8, units="cm",res=155,pointsize = 20)
autoplot(pca_outliers, data = strata_5pops, colour = 'pop'  ) + theme_bw() + scale_color_manual(values=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()




# vcf  neutral snps ####

system2("vcftools", 
        args = c("--vcf", "./output/batch_1_JTR_JPA.vcf",
                 "--stdout", "--recode",
                 "--recode-INFO-all", "--exclude", "./output/outliers_PCAdapt_Bayescan.txt",
                 ">", "./output/neutral_snps_JTR_JPA.vcf"))


vcf_JTR_JPA_neutral = read.vcf("./output/neutral_snps_JTR_JPA.vcf", which.loci = 1:1e6)
dim(vcf_JTR_JPA_neutral) #  98  indv, 17256 SNPs (18879 snps - 1623 outlier )


genind_JTR_JPA_neutral <- loci2genind(vcf_JTR_JPA_neutral)
summary_JTR_JPA_neutral = summary(genind_JTR_JPA_neutral)

str(summary_JTR_JPA_neutral)
summary_JTR_JPA_neutral$n.by.pop
mean(summary_JTR_JPA_neutral$Hobs) # 0.1746265
mean(summary_JTR_JPA_neutral$Hexp) # 0.2262396
min(summary_JTR_JPA_neutral$loc.n.all)
max(summary_JTR_JPA_neutral$loc.n.all)
mean(summary_JTR_JPA_neutral$NA.perc) # 30.12102


# strata 5 pops #### 

JTR_JPA_strata_5pops <- read.table("./input/strata_JTR_JPA_5pops.txt", header=T)
head(JTR_JPA_strata_5pops)
dim(JTR_JPA_strata_5pops)

strata(genind_JTR_JPA_neutral) <- JTR_JPA_strata_5pops
genind_JTR_JPA_neutral

library("mmod")
setPop(genind_JTR_JPA_neutral) = ~pop
str(genind_JTR_JPA_neutral)

summary_JTR_JPA_5pops = summary(genind_JTR_JPA_neutral) # 98 indv, 5 pops
summary_JTR_JPA_5pops$n.by.pop


# Df genind ordered ####

df_JTR_JPA_neutral = genind2df(genind_JTR_JPA_neutral, sep = "/")
df_JTR_JPA_neutral[1:15,1:5]

write.csv(df_JTR_JPA_neutral, file = "./output/df_JTR_JPA_neutral.csv", row.names = T)
JTR_JPA_neutral_ordered <- read.table("./output/df_JTR_JPA_neutral_ordered.csv", header=T, sep = ",")
JTR_JPA_neutral_ordered[1:5,1:5]
dim(JTR_JPA_neutral_ordered)

JTR_JPA_neutral_genind_ordered = df2genind(JTR_JPA_neutral_ordered, sep = "/")

JTR_JPA_strata_5pops_ordered <- read.table("./input/strata_JTR_JPA_5pops_ordered.txt", header=T)
head(JTR_JPA_strata_5pops_ordered)
dim(JTR_JPA_strata_5pops_ordered)

strata(JTR_JPA_neutral_genind_ordered) <- JTR_JPA_strata_5pops_ordered
JTR_JPA_neutral_genind_ordered

library("mmod")
setPop(JTR_JPA_neutral_genind_ordered) = ~pop
str(JTR_JPA_neutral_genind_ordered)

summary_JTR_JPA_5pops_ordered = summary(JTR_JPA_neutral_genind_ordered) # 98 indv, 5 pops
summary_JTR_JPA_5pops_ordered$n.by.pop



#  Ar neutral ####

library(hierfstat)

JTR_JPA_hierfstat_neutral = genind2hierfstat(genind_JTR_JPA_neutral)
dim(JTR_JPA_hierfstat_neutral)
JTR_JPA_hierfstat_neutral[1:5,1:5]

allelic.richness_JTR_JPA_neutral = allelic.richness(JTR_JPA_hierfstat_neutral)
str(allelic.richness_JTR_JPA_neutral)
dim(allelic.richness_JTR_JPA_neutral$Ar)

AR_matrix_JTR_JPA_neutral = as.data.frame(allelic.richness_JTR_JPA_neutral$Ar)
str(AR_matrix_JTR_JPA_neutral)
dim(AR_matrix_JTR_JPA_neutral)
write.csv(AR_matrix_JTR_JPA_neutral, file = "./output/AR_matrix_JTR_JPA_neutral.csv")


# 2 pops #

JTR_JPA_strata <- read.table("./input/strata_JTR_JPA.txt", header=T)
head(JTR_JPA_strata)
dim(JTR_JPA_strata)

strata(genind_JTR_JPA_neutral) <- JTR_JPA_strata
genind_JTR_JPA_neutral

library("mmod")
setPop(genind_JTR_JPA_neutral) = ~pop
str(genind_JTR_JPA_neutral)


summary_JTR_JPA_neutral = summary(genind_JTR_JPA_neutral) # 98 indv, 2 pops
str(summary_JTR_JPA_neutral)

JTR_JPA_hierfstat_2pops_neutral = genind2hierfstat(genind_JTR_JPA_neutral)
dim(JTR_JPA_hierfstat_2pops_neutral)
JTR_JPA_hierfstat_2pops_neutral[1:5,1:5]

allelic.richness_JTR_JPA_neutral_2pops = allelic.richness(JTR_JPA_hierfstat_2pops_neutral)
str(allelic.richness_JTR_JPA_neutral_2pops)
dim(allelic.richness_JTR_JPA_neutral_2pops$Ar)

AR_matrix_JTR_JPA_neutral_2pops = as.data.frame(allelic.richness_JTR_JPA_neutral_2pops$Ar)
str(AR_matrix_JTR_JPA_neutral_2pops)
dim(AR_matrix_JTR_JPA_neutral_2pops)
write.csv(AR_matrix_JTR_JPA_neutral_2pops, file = "./output/AR_matrix_JTR_JPA_neutral_2pops.csv")




# inbreeding neutral ####
# adegenet #


JTR_genind = seppop(genind_JTR_JPA_neutral)$JTR
Fis_JTR = inbreeding(JTR_genind)
Fis_mean_JTR = sapply(Fis_JTR, mean)
mean(Fis_mean_JTR)

JPA_genind = seppop(genind_JTR_JPA_neutral)$JPA
Fis_JPA = inbreeding(JPA_genind)
Fis_mean_JPA = sapply(Fis_JPA, mean)
mean(Fis_mean_JPA)




# HO HE neutral ####

# adegenet #
genind_JTR_JPA_neutral$pop
# Levels: JPA_SP JPA_AM JTR_TC JTR_GI JTR_NI

JPA_SP_genind_neutral = seppop(genind_JTR_JPA_neutral)$JPA_SP
summary_JPA_SP_genind_neutral = summary(JPA_SP_genind_neutral)
str(summary_JPA_SP_genind_neutral)
mean(summary_JPA_SP_genind_neutral$Hobs, na.rm=T)
mean(summary_JPA_SP_genind_neutral$Hexp, na.rm=T)

JPA_AM_genind_neutral = seppop(genind_JTR_JPA_neutral)$JPA_AM
summary_JPA_AM_genind_neutral = summary(JPA_AM_genind_neutral)
str(summary_JPA_AM_genind_neutral)
mean(summary_JPA_AM_genind_neutral$Hobs, na.rm=T)
mean(summary_JPA_AM_genind_neutral$Hexp, na.rm=T)

JTR_TC_genind_neutral = seppop(genind_JTR_JPA_neutral)$JTR_TC
summary_JTR_TC_genind_neutral = summary(JTR_TC_genind_neutral)
str(summary_JTR_TC_genind_neutral)
mean(summary_JTR_TC_genind_neutral$Hobs, na.rm=T)
mean(summary_JTR_TC_genind_neutral$Hexp, na.rm=T)

JTR_GI_genind_neutral = seppop(genind_JTR_JPA_neutral)$JTR_GI
summary_JTR_GI_genind_neutral = summary(JTR_GI_genind_neutral)
str(summary_JTR_GI_genind_neutral)
mean(summary_JTR_GI_genind_neutral$Hobs, na.rm=T)
mean(summary_JTR_GI_genind_neutral$Hexp, na.rm=T)

JTR_NI_genind_neutral = seppop(genind_JTR_JPA_neutral)$JTR_NI
summary_JTR_NI_genind_neutral = summary(JTR_NI_genind_neutral)
str(summary_JTR_NI_genind_neutral)
mean(summary_JTR_NI_genind_neutral$Hobs, na.rm=T)
mean(summary_JTR_NI_genind_neutral$Hexp, na.rm=T)


# 2 pops #

genind_JTR_JPA_neutral$pop
# Levels: JPA JTR

JPA_genind_neutral = seppop(genind_JTR_JPA_neutral)$JPA
summary_JPA_genind_neutral = summary(JPA_genind_neutral)
mean(summary_JPA_genind_neutral$Hobs, na.rm=T)
mean(summary_JPA_genind_neutral$Hexp, na.rm=T)

JTR_genind_neutral = seppop(genind_JTR_JPA_neutral)$JTR
summary_JTR_genind_neutral = summary(JTR_genind_neutral)
mean(summary_JTR_genind_neutral$Hobs, na.rm=T)
mean(summary_JTR_genind_neutral$Hexp, na.rm=T)




# DAPC ####

grp_JTR_JPA_neutral <- find.clusters(genind_JTR_JPA_neutral, max.n.clust=20) # retain 100 PCs, 2 clusters
dapc_JTR_JPA_neutral <- dapc(genind_JTR_JPA_neutral, grp_JTR_JPA_neutral$grp_JTR_JPA_neutral)
# 98/3= 33 PCs,  $var = 0.4315677, a-score = 
# retain 1 PC (optimal),  $var = 0.01901853,  a-score = 
# 2 PCs, $var = 0.03510191,  a-score = 

head(dapc_JTR_JPA_neutral) # 


# a-score  ##

temp <- a.score(dapc_JTR_JPA_neutral)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score
temp$mean # 


png("./figs/a-score_optimisations_JTR_JPA_neutral.png", width=50, height=30, units="cm",res=155,pointsize = 36)
temp <- optim.a.score(dapc_JTR_JPA_neutral) # 1 PC optimal
dev.off()

myCol5 <- c("darkblue", "lightblue", "green", "darkgreen", "lightgreen") 

png("./figs/PCA_JTR_JPA_neutral_2PCs.png", width=30, height=20, units="cm",res=155,pointsize = 24)
scatter(dapc_JTR_JPA_neutral, posi.da="bottomright", bg="white", pch=17:22)
dev.off()

png("./figs/DAPC_JTR_JPA_neutral_2PCs.png", width=30, height=20, units="cm",res=155,pointsize = 18)
scatter(dapc_JTR_JPA_neutral, scree.da=FALSE, bg="white", pch=20, cell=1.5, cstar=0, 
        cex=3,clab=0, leg=T, txt.leg=paste(c("JPA_SP", "JPA_AM", "JTR_TC-II", "JTR_GI", "JTR_NI")), 
        posi.leg = "bottomleft", col=myCol5)
dev.off()


#  memb prob ####


png("./figs/memb_prob_JTR_JPA_neutral_2Pcs_spp.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_neutral, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()

png("./figs/memb_prob_JTR_JPA_neutral_1Pcs_spp.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_neutral, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()

png("./figs/memb_prob_JTR_JPA_neutral_1Pcs_2pops.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_neutral, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00", "#0d47a1"))
dev.off()

# select admixed ind ##

admix <- which(apply(dapc_JTR_JPA_neutral$posterior,1, function(e) all(e<0.999)))
admix

png("./figs/memb_prob_JTR_JPA_neutral_1Pcs_admixed.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_neutral,subset=admix, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#0d47a1"))
dev.off()



# DAPC ordered ####

grp_JTR_JPA_neutral_ordered <- find.clusters(JTR_JPA_neutral_genind_ordered, max.n.clust=20) # retain 100 PCs, 2 clusters
dapc_JTR_JPA_neutral_ordered <- dapc(JTR_JPA_neutral_genind_ordered, grp_JTR_JPA_neutral_ordered$grp_JTR_JPA_neutral_ordered)
# retain 1 PC (optimal),  $var = 0.01901853,  a-score = 0.3090909

head(dapc_JTR_JPA_neutral_ordered) # 

# a-score  ##

temp <- a.score(dapc_JTR_JPA_neutral_ordered)
names(temp)
temp$tab[1:2,1:2]
temp$pop.score
temp$mean # 


png("./figs/PCA_JTR_JPA_neutral_2PCs.png", width=30, height=20, units="cm",res=155,pointsize = 24)
scatter(dapc_JTR_JPA_neutral_ordered, posi.da="bottomright", bg="white", pch=17:22)
dev.off()


#  memb prob ordered ####


png("./figs/memb_prob_JTR_JPA_neutral_1Pc_ordered.png", width=50, height=30, units="cm",res=155,pointsize = 36)
compoplot(dapc_JTR_JPA_neutral_ordered, lab="", posi=list(x=12,y=-.01), cleg=.7, include.origin = TRUE, col=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()




# PCA neutral ####

library(hierfstat)

JTR_JPA_neutral_hierfstat = genind2hierfstat(genind_JTR_JPA_neutral)
dim(JTR_JPA_neutral_hierfstat)
JTR_JPA_neutral_hierfstat[1:5,1:5]

write.csv(JTR_JPA_neutral_hierfstat, file = "./output/JTR_JPA_neutral_hierfstat.csv", row.names = F)
# replace "NAs" for "0"
JTR_JPA_neutral_hierfstat_PCA = read.csv("./output/JTR_JPA_neutral_hierfstat_PCA.csv", header = T)
dim(JTR_JPA_neutral_hierfstat_PCA) # 98 17257

JTR_JPA_neutral_hierfstat_PCA = JTR_JPA_neutral_hierfstat_PCA[,2:17257]
JTR_JPA_neutral_hierfstat_PCA[1:5,1:5]


JTR_JPA_indv_list = read.csv("./input/JTR_JPA_indv_list.txt", header = F)
JTR_JPA_indv_list = as.matrix(JTR_JPA_indv_list)
head(JTR_JPA_indv_list) 
JTR_JPA_neutral_hierfstat_PCA$id = JTR_JPA_indv_list
dim(JTR_JPA_neutral_hierfstat_PCA)
JTR_JPA_neutral_hierfstat_PCA[1:5,17250:17257] 
JTR_JPA_neutral_hierfstat_PCA = JTR_JPA_neutral_hierfstat_PCA[,c(17257,1:17256)]
str(JTR_JPA_neutral_hierfstat_PCA)
JTR_JPA_neutral_hierfstat_PCA[1:5,1:5]

write.csv(JTR_JPA_neutral_hierfstat_PCA, file = "./output/JTR_JPA_neutral_hierfstat_PCA.csv", row.names = F)
JTR_JPA_neutral_hierfstat_PCA = read.csv("./output/JTR_JPA_neutral_hierfstat_PCA.csv", header = T, row.names = 1)
JTR_JPA_neutral_hierfstat_PCA[1:5,1:5]


# row.names(JTR_JPA_outliers_hierfstat_PCA) = JTR_JPA_outliers_hierfstat_PCA[,1]

pca = prcomp(JTR_JPA_neutral_hierfstat_PCA)
summary(pca)
# PC1
# Standard deviation     893.1936
# Proportion of Variance   0.2988
# PC2
# Standard deviation    368.3031
# Proportion of Variance  0.0508



png("./figs/pca_neutral_JTR_JPA.png", width=40, height=30, units="cm",res=155,pointsize = 26)
plot(pca$x, pch=20, col="blue") #  xlim=c(-11,9), ylim=c(-4.5,4.5)
text(pca$x, rownames(pca$x), pos=3, cex=0.6)
dev.off()


strata_5pops = read.csv("./input/strata_JTR_JPA_5pops.txt")
library(ggfortify)
png("./figs/pca_neutral_JTR_JPA_5pops.png", width=15, height=8, units="cm",res=155,pointsize = 20)
autoplot(pca, data = strata_5pops, colour = 'pop'  ) + theme_bw() + scale_color_manual(values=c("#ef6c00","#ffb74d", "#0d47a1","#1e88e5","#90caf9"))
dev.off()



# downsample SNPs ####

# vcf  JTR_JPA batch 1 ####
vcf_JTR_JPA = read.vcf("./output/batch_1_JTR_JPA.vcf", which.loci = 1:1e6)
dim(vcf_JTR_JPA) #  98  indv, 18879 SNPs

# randomly select 100 SNPs

# cd /Users/jc471612/Documents/James_Cook/Lobsters/3_Jasus/7_JTR_JPA/R/output
# gshuf -n 100 ./batch_1_JTR_JPA.vcf > ./batch_1_JTR_JPA_downsampled_1.vcf


# select neutral SNPs only for migrate???




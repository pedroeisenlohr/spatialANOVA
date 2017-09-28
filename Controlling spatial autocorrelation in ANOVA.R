##########################################################################################################
### Variance partitioning as a tool for controlling the type I error in ANOVA routine: a suggested R code.
##########################################################################################################
### After Eisenlohr (2014)
### Released on: 10 August, 2017

### If you identify any problem or have any suggestion to improve this code, please contact the author of the manuscript (pedrov.eisenlohr@gmail.com).
### This code was prepared and edited by Danilo R.M. Neves (University of Leeds) and Pedro V. Eisenlohr (Universidade Federal de Minas Gerais), with the help of Mário José Marques Azevedo (Universidade Estadual de Campinas), considering the following references:
# 1. Legendre, P., D. Borcard and D. W. Roberts. 2012. Variation partitioning involving orthogonal spatial eigenfunction submodels. Ecology 93: 1234-1240.
# 2. Borcard, D., F. Gillet & P. Legendre. 2011. Numerical Ecology with R. Springer, New York, 302p.
# 3. Peres-Neto, P. R. and P. Legendre. 2010. Estimating and controlling for spatial structure in the study of ecological communities. Global Ecology and Biogeography 19: 174-184.
# 4. Blanchet F. G., P. Legendre, and D. Borcard. 2008. Forward selection of explanatory variables. Ecology 89: 2623-2632.
# 5. Dray, S., P. Legendre and P. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbor matrices (PCNM). Ecological Modelling 196: 483-493.
# 6. Legendre, P. and E. Gallagher. 2001. Ecologically meaningful transformations for ordination of species data. Oecologia 129: 271-280.
# 7. Anderson, M. J. & P. Legendre. 1999. An empirical comparison of permutation methods for tests of partial regression coefficients in a linear model. Journal of Statistical Computation and Simulation 62: 271-303.

### Preparing data matrices:
# X matrix: full model (without latitude and longitude). VarDep is the response variable
# ll matrix: latitude and longitude
# Y matrix: only response variable (s) (VarDep or, in case of MANOVA, VarDep1, VarDep 2 etc.)
# Examples of such matrices can be found assessing Supplementary Material of Eisenlohr (2014).

# Define the working directory.
# Each user should adjust this!
setwd("C:/Users/pedro/Dropbox/Ferramentas LabEc/Matrizes")
getwd()


# If you do not have the packages installed, please use the following commands:
install.packages("ape") 
install.packages("packfor", repos="http://R-Forge.R-project.org")
install.packages("spacemakeR", repos="http://R-Forge.R-project.org")
install.packages("ade4")
install.packages("spdep") 
install.packages("vegan")
install.packages("tripack")
install.packages("AEM", repos="http://R-Forge.R-project.org")
install.packages("PCNM", repos="http://R-Forge.R-project.org")
# PCNM and/or AEM package have/has experienced problems with this link. If it happens, please 
# install it manually (.zip files) from another source.

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ape)
library(packfor)
library(spacemakeR)
library(ade4)
library(spdep)
library(vegan)
library(tripack)
library(AEM)
library(PCNM)

### Importing and verifying data matrices:
X<-read.table(file.choose(),row.names=1,header=T,sep=",") 
dim (X)
edit(X)
ll<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim (ll)
edit(ll)
Y<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim (Y)
edit(Y)

### Generating and testing the significance of spatial variables 
### (Moran's Eigenvector Maps - MEMs), and choosing the best MEM Type:

################ MEM TYPE 1: Classic PCNM #############################
## Construct the PCNM variables
xy.d1 <- dist(ll)
spp.PCNM.auto <- PCNM(xy.d1)
# Truncation distance
(dmin <- spp.PCNM.auto$thresh)
# Number of eigenvalues
(nb.ev <- length(spp.PCNM.auto$values))
# Expected value of I, no spatial correlation
spp.PCNM.auto$expected_Moran
# Moran's I of the PCNM variables
spp.PCNM.auto$Moran_I
# Eigenfunctions with positive spatial correlation
(select <- which(spp.PCNM.auto$Moran_I$Positive == TRUE))
# Number of PCNM with I > E(I)
length(select)
spp.PCNM.pos <- as.data.frame(spp.PCNM.auto$vectors)[,select]
write.table(spp.PCNM.pos,"PCNMpresel.csv")

### Checking for the significance of this MEM subset:
lm1<-lm(VarDep~spp.PCNM.auto$vectors[,select], data=Y)
summary(lm1)

### According to Blanchet et al. (2008): "If, and only if, the global test is
### significant, one can proceed with forward selection"

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm1<-RsquareAdj(lm1)$adj.r.squared
r2.lm1
sel.lm1<-forward.sel(Y,spp.PCNM.pos,adjR2thresh=r2.lm1)
sel.lm1
# Write the significant PCNMs to a new object
(nb.sig.PCNM <- nrow(sel.lm1)) #Number of signif. PCNM
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(sel.lm1[,2]))
# Write the significant PCNMs to a new object
PCNM.red <- spp.PCNM.pos[,c(PCNM.sign)]
write.table(PCNM.red,"PCNMsel.csv")
### Here it's important to register the R2a after PCNM selection. This value will
### be compared with the R2a obtained for alternative spatial models (MEMs 2-16;
### see below), which should also be recorded.

################ MEM TYPE 2 #############################
# Delaunay triangulation. No weighting matrix (binary weights).
ll.temp2 <- tri2nb(ll)
ll.transf2 = nb2listw(ll.temp2)
mem2<-scores.listw(ll.transf2) 
colnames(mem2$vectors)<-paste("MEMT2",1:ncol(mem2$vectors))
colnames(mem2$vectors)
mem.teste2<-test.scores(mem2,ll.transf2,1000)
mem.teste2
sel.mem2<-mem2
sel.mem2$vectors<-sel.mem2$vectors[,mem.teste2[,2]<0.05]
dim(sel.mem2$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos2<-which(mem.teste2[,1]>-1/(nrow(mem2$vectors)-1))
MEM.pos2<-mem2$vectors[,MEM.Moran.pos2]
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig2<-MEM.Moran.pos2[which(mem.teste2[MEM.Moran.pos2,2]<=0.05)]
MEM.pos.sig2<-mem2$vectors[,MEM.Moran.pos.sig2]

### Checking for the significance of this MEM subset:
lm2<-lm(VarDep~MEM.pos.sig2, data=Y)
summary(lm2)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm2<-RsquareAdj(lm2)$adj.r.squared
r2.lm2
sel.lm2<-forward.sel(Y,MEM.pos.sig2,adjR2thresh=r2.lm2)
sel.lm2
espaciais<-colnames(mem2$vectors)%in%sel.lm2$variables
sel.values<-sel.mem2$values[sel.lm2$order]
sel.vectors2<-sel.mem2$vectors[,sel.lm2$order]
write.table(sel.vectors2,"final.sel_MEM-T2.csv")

################ MEM TYPE 3 #############################
# Delaunay triangulation. "Minmax" matrix. 
### The ?minmax? style is based on Kelejian and Prucha (2010),
### and divides the weights by the minimum of the maximum row sums
### and maximum column sums of the input weights.
ll.temp3<-tri2nb(ll) 
ll.transf3 = nb2listw(ll.temp3, style = "minmax")
mem3<-scores.listw(ll.transf3) 
colnames(mem3$vectors)<-paste("MEMT3",1:ncol(mem3$vectors))
colnames(mem3$vectors)
mem.teste3<-test.scores(mem3,ll.transf3,1000)
mem.teste3
sel.mem3<-mem3
sel.mem3$vectors<-sel.mem3$vectors[,mem.teste3[,2]<0.05]
dim(sel.mem3$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos3<-which(mem.teste3[,1]>-1/(nrow(mem3$vectors)-1))
MEM.pos3<-mem3$vectors[,MEM.Moran.pos3]
dim(MEM.pos3)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig3<-MEM.Moran.pos3[which(mem.teste3[MEM.Moran.pos3,2]<=0.05)]
MEM.pos.sig3<-mem3$vectors[,MEM.Moran.pos.sig3]
dim(MEM.pos.sig3)

### Checking for the significance of this MEMs subset:
lm3<-lm(VarDep~MEM.pos.sig3, data=Y)
summary(lm3)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm3<-RsquareAdj(lm3)$adj.r.squared
r2.lm3
sel.lm3<-forward.sel(Y,MEM.pos.sig3,adjR2thresh=r2.lm3)
sel.lm3
espaciais<-colnames(mem3$vectors)%in%sel.lm3$variables
sel.values<-sel.mem3$values[sel.lm3$order]
sel.vectors3<-sel.mem3$vectors[,sel.lm3$order]
write.table(sel.vectors3,"final.sel_MEM-T3.csv")

################ MEM TYPE 4 #############################
# Delaunay triangulation. "B" matrix.
ll.temp4<-tri2nb(ll) 
ll.transf4 = nb2listw(ll.temp4, style = "B")
mem4<-scores.listw(ll.transf4) 
colnames(mem4$vectors)<-paste("MEMT4",1:ncol(mem4$vectors))
colnames(mem4$vectors)
mem.teste4<-test.scores(mem4,ll.transf4,1000)
sel.mem4<-mem4
sel.mem4$vectors<-sel.mem4$vectors[,mem.teste4[,2]<0.05]
dim(sel.mem4$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos4<-which(mem.teste4[,1]>-1/(nrow(mem4$vectors)-1))
MEM.pos4<-mem4$vectors[,MEM.Moran.pos4]
dim(MEM.pos4)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig4<-MEM.Moran.pos4[which(mem.teste4[MEM.Moran.pos4,2]<=0.05)]
MEM.pos.sig4<-mem4$vectors[,MEM.Moran.pos.sig4]
dim(MEM.pos.sig4)

### Checking for the significance of this MEMs subset:
lm4<-lm(VarDep~MEM.pos.sig4, data=Y)
summary(lm4)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm4<-RsquareAdj(lm4)$adj.r.squared
r2.lm4
sel.lm4<-forward.sel(Y,MEM.pos.sig4,adjR2thresh=r2.lm4)
sel.lm4
espaciais<-colnames(mem4$vectors)%in%sel.lm4$variables
sel.values<-sel.mem4$values[sel.lm4$order]
sel.vectors4<-sel.mem4$vectors[,sel.lm4$order]
write.table(sel.vectors4,"final.sel_MEM-T4.csv")

################ MEM TYPE 5 #############################
# Delaunay triangulation. "C" matrix.
ll.temp5<-tri2nb(ll) 
ll.transf5 = nb2listw(ll.temp5, style = "C")
mem5<-scores.listw(ll.transf5) 
colnames(mem5$vectors)<-paste("MEMT5",1:ncol(mem5$vectors))
colnames(mem5$vectors)
mem.teste5<-test.scores(mem5,ll.transf5,1000)
sel.mem5<-mem5
sel.mem5$vectors<-sel.mem5$vectors[,mem.teste5[,2]<0.05]
dim(sel.mem5$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos5<-which(mem.teste5[,1]>-1/(nrow(mem5$vectors)-1))
MEM.pos5<-mem5$vectors[,MEM.Moran.pos5]
dim(MEM.pos5)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig5<-MEM.Moran.pos5[which(mem.teste5[MEM.Moran.pos5,2]<=0.05)]
MEM.pos.sig5<-mem5$vectors[,MEM.Moran.pos.sig5]
dim(MEM.pos.sig5)

### Checking for the significance of this MEMs subset:
lm5<-lm(VarDep~MEM.pos.sig5, data=Y)
summary(lm5)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm5<-RsquareAdj(lm5)$adj.r.squared
r2.lm5
sel.lm5<-forward.sel(Y,MEM.pos.sig5,adjR2thresh=r2.lm5)
sel.lm5
espaciais<-colnames(mem5$vectors)%in%sel.lm5$variables
sel.values<-sel.mem5$values[sel.lm5$order]
sel.vectors5<-sel.mem5$vectors[,sel.lm5$order]
write.table(sel.vectors5,"final.sel_MEM-T5.csv")

################ MEM TYPE 6 #############################
# Delaunay triangulation. "U" matrix.
ll.temp6<-tri2nb(ll) 
ll.transf6 = nb2listw(ll.temp6, style = "U")
mem6<-scores.listw(ll.transf6) 
colnames(mem6$vectors)<-paste("MEMT6",1:ncol(mem6$vectors))
colnames(mem6$vectors)
mem.teste6<-test.scores(mem6,ll.transf6,1000)
sel.mem6<-mem6
sel.mem6$vectors<-sel.mem6$vectors[,mem.teste6[,2]<0.05]
dim(sel.mem6$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos6<-which(mem.teste6[,1]>-1/(nrow(mem6$vectors)-1))
MEM.pos6<-mem6$vectors[,MEM.Moran.pos6]
dim(MEM.pos6)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig6<-MEM.Moran.pos6[which(mem.teste6[MEM.Moran.pos6,2]<=0.05)]
MEM.pos.sig6<-mem6$vectors[,MEM.Moran.pos.sig6]
dim(MEM.pos.sig6)

### Checking for the significance of this MEMs subset:
lm6<-lm(VarDep~MEM.pos.sig6, data=Y)
summary(lm6)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm6<-RsquareAdj(lm6)$adj.r.squared
r2.lm6
sel.lm6<-forward.sel(Y,MEM.pos.sig6,adjR2thresh=r2.lm6)
sel.lm6
espaciais<-colnames(mem6$vectors)%in%sel.lm6$variables
sel.values<-sel.mem6$values[sel.lm6$order]
sel.vectors6<-sel.mem6$vectors[,sel.lm6$order]
write.table(sel.vectors6,"final.sel_MEM-T6.csv")

################ MEM TYPE 7 #############################
# Delaunay triangulation. "S" matrix.
ll.temp7<-tri2nb(ll) 
ll.transf7 = nb2listw(ll.temp7, style = "S")
mem7<-scores.listw(ll.transf7) 
colnames(mem7$vectors)<-paste("MEMT7",1:ncol(mem7$vectors))
colnames(mem7$vectors)
mem.teste7<-test.scores(mem7,ll.transf7,1000)
sel.mem7<-mem7
sel.mem7$vectors<-sel.mem7$vectors[,mem.teste7[,2]<0.05]
dim(sel.mem7$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos7<-which(mem.teste7[,1]>-1/(nrow(mem7$vectors)-1))
MEM.pos7<-mem7$vectors[,MEM.Moran.pos7]
dim(MEM.pos7)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig7<-MEM.Moran.pos7[which(mem.teste7[MEM.Moran.pos7,2]<=0.05)]
MEM.pos.sig7<-mem7$vectors[,MEM.Moran.pos.sig7]
dim(MEM.pos.sig7)

### Checking for the significance of this MEMs subset:
lm7<-lm(VarDep~MEM.pos.sig7, data=Y)
summary(lm7)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm7<-RsquareAdj(lm7)$adj.r.squared
r2.lm7
sel.lm7<-forward.sel(Y,MEM.pos.sig7,adjR2thresh=r2.lm7)
sel.lm7
espaciais<-colnames(mem7$vectors)%in%sel.lm7$variables
sel.values<-sel.mem7$values[sel.lm7$order]
sel.vectors7<-sel.mem7$vectors[,sel.lm7$order]
write.table(sel.vectors7,"final.sel_MEM-T7.csv")

################ MEM TYPE 8 #############################
# Gabriel graph
ll.temp8 <- graph2nb(gabrielneigh(as.matrix(ll)), sym=TRUE)
ll.transf8 = nb2listw(ll.temp8)
mem8<-scores.listw(ll.transf8) 
colnames(mem8$vectors)<-paste("MEMT8",1:ncol(mem8$vectors))
colnames(mem8$vectors)
mem.teste8<-test.scores(mem8,ll.transf8,1000)
sel.mem8<-mem8
sel.mem8$vectors<-sel.mem8$vectors[,mem.teste8[,2]<0.05]
dim(sel.mem8$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos8<-which(mem.teste8[,1]>-1/(nrow(mem8$vectors)-1))
MEM.pos8<-mem8$vectors[,MEM.Moran.pos8]
dim(MEM.pos8)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig8<-MEM.Moran.pos8[which(mem.teste8[MEM.Moran.pos8,2]<=0.05)]
MEM.pos.sig8<-mem8$vectors[,MEM.Moran.pos.sig8]
dim(MEM.pos.sig8)

### Checking for the significance of this MEMs subset:
lm8<-lm(VarDep~MEM.pos.sig8, data=Y)
summary(lm8)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm8<-RsquareAdj(lm8)$adj.r.squared
r2.lm8
sel.lm8<-forward.sel(Y,MEM.pos.sig8,adjR2thresh=r2.lm8)
sel.lm8
espaciais<-colnames(mem8$vectors)%in%sel.lm8$variables
sel.values<-sel.mem8$values[sel.lm8$order]
sel.vectors8<-sel.mem8$vectors[,sel.lm8$order]
write.table(sel.vectors8,"final.sel_MEM-T8.csv")

################ MEM TYPE 9 #############################
# Relative neighbourhood
ll.temp9 <- graph2nb(relativeneigh(as.matrix(ll)), sym=TRUE)
ll.transf9 = nb2listw(ll.temp9)
mem9<-scores.listw(ll.transf9) 
colnames(mem9$vectors)<-paste("MEMT9",1:ncol(mem9$vectors))
colnames(mem9$vectors)
mem.teste9<-test.scores(mem9,ll.transf9,1000)
sel.mem9<-mem9
sel.mem9$vectors<-sel.mem9$vectors[,mem.teste9[,2]<0.05]
dim(sel.mem9$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos9<-which(mem.teste9[,1]>-1/(nrow(mem9$vectors)-1))
MEM.pos9<-mem9$vectors[,MEM.Moran.pos9]
dim(MEM.pos9)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig9<-MEM.Moran.pos9[which(mem.teste9[MEM.Moran.pos9,2]<=0.05)]
MEM.pos.sig9<-mem9$vectors[,MEM.Moran.pos.sig9]
dim(MEM.pos.sig9)

### Checking for the significance of this MEMs subset:
lm9<-lm(VarDep~MEM.pos.sig9, data=Y)
summary(lm9)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm9<-RsquareAdj(lm9)$adj.r.squared
r2.lm9
sel.lm9<-forward.sel(Y,MEM.pos.sig9,adjR2thresh=r2.lm9)
sel.lm9
espaciais<-colnames(mem9$vectors)%in%sel.lm9$variables
sel.values<-sel.mem9$values[sel.lm9$order]
sel.vectors9<-sel.mem9$vectors[,sel.lm9$order]
write.table(sel.vectors9,"final.sel_MEM-T9.csv")

################ MEM TYPE 10 #############################
# Minimum spanning tree
ll.temp10 <- mst.nb(dist(ll))
ll.transf10 = nb2listw(ll.temp10)
mem10<-scores.listw(ll.transf10) 
colnames(mem10$vectors)<-paste("MEMT10",1:ncol(mem10$vectors))
colnames(mem10$vectors)
mem.teste10<-test.scores(mem10,ll.transf10,1000)
sel.mem10<-mem10
sel.mem10$vectors<-sel.mem10$vectors[,mem.teste10[,2]<0.05]
dim(sel.mem10$vectors)
#MEM with positive spatial autocorrelation
MEM.Moran.pos10<-which(mem.teste10[,1]>-1/(nrow(mem10$vectors)-1))
MEM.pos10<-mem10$vectors[,MEM.Moran.pos10]
dim(MEM.pos10)
#MEM with positive and significant spatial autocorrelation
MEM.Moran.pos.sig10<-MEM.Moran.pos10[which(mem.teste10[MEM.Moran.pos10,2]<=0.05)]
MEM.pos.sig10<-mem10$vectors[,MEM.Moran.pos.sig10]
dim(MEM.pos.sig10)

### Checking for the significance of this MEMs subset:
lm10<-lm(VarDep~MEM.pos.sig10, data=Y)
summary(lm10)

## Since the analysis is significant, compute the adjusted R2 and run a forward
## selection of the PCNM variables.
r2.lm10<-RsquareAdj(lm10)$adj.r.squared
r2.lm10
sel.lm10<-forward.sel(Y,MEM.pos.sig10,adjR2thresh=r2.lm10)
sel.lm10
espaciais<-colnames(mem10$vectors)%in%sel.lm10$variables
sel.values<-sel.mem10$values[sel.lm10$order]
sel.vectors10<-sel.mem10$vectors[,sel.lm10$order]
write.table(sel.vectors10,"final.sel_MEM-T10.csv")

############################################################################
### The following commands should be used only for checking the ANOVA/MANOVA
### results when spatial structures are not considered.#####################
############################################################################
lm<-lm(VarDep~treat, data=X)
anova(lm)

### Preparing selected variables for variance partitioning:
### Please change * by the number of the best spatial model.
### In case of the classical PCNM was the best MEM type, you should skip this step.
spatial<-colnames(mem*$vectors)%in%sel.lm*$variables
sel.values<-sel.mem*$values[sel.lm*$order]
sel.vectors<-sel.mem*$vectors[,sel.lm*$order]
dim(sel.vectors)
write.table(sel.vectors,"mem_selected.csv")

treat <- model.matrix(VarDep~treat, data=X)[ ,-1]

### Variance partitioning:
### In case of classical PCNM (MEM Type 1) elected as the best MEM type, please change 
### sel.vectors by PCNM.red from now on.
all.varpart<-varpart(Y,treat,sel.vectors) 
all.varpart
plot(all.varpart)

### Testing the significance of the treatment (fraction [a]), after considering the effect of selected MEMs:
amb.rda<-rda(Y,treat,sel.vectors)
teste.amb<-anova(amb.rda)
teste.amb

### In case of significance, you should perform post hoc tests (Tukey etc.) using the selected MEMs as covariables.

#################################################################################
### If you would like to clean all the steps above, please run:
### PLEASE NOTE: THE FOLLOWING COMMAND WILL ERASE ALL THE MEMORY OF THIS ROUTINE!
rm(list=ls(all=TRUE))
#################################################################################


##########################################################################################
######################################### TO DO LIST (Pedro's list): ####################################
1) Make some adjustments to PCNM routine.
2) Fix some errors on MEM11-MEM16 (which were removed from this file due to such errors).
3) Change 'lm' by 'lmp' in global tests.
4) Implement post hoc test.
##########################################################################################
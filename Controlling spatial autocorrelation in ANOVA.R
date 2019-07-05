##########################################################################################################
### ANOVA controlling spatial structures ###
##########################################################################################################
### After Eisenlohr (2014) https://doi.org/10.1007/s40415-014-0064-3
### Released on: 04 July, 2019


# Define the working directory.
# Each user should adjust this!
setwd(choose.dir())
getwd()

# Load the required packages
library(adespatial)

### Importing and verifying data matrices:
matrix<-read.table(file.choose(),row.names=1,header=T,sep=",") #Include treat (groups) and Y (response variable)
dim (matrix)
View(matrix)
ll<-read.table(file.choose(),row.names=1,header=T,sep=",") #Longitude and Latitude
dim (ll)
View(ll)

Y <- as.data.frame(matrix$Y)
View(Y)

### Generating and testing the significance of spatial variables 
### (Moran's Eigenvector Maps - MEMs), and choosing the best MEM Type:

##### OPTIMIZING THE SELECTION OF SMW #####
### The function listw.candidates is used to build the spatial weighting matrices that
### we want to test and compare (with the listw.select function).
### I strongly recommend a careful reading on Bauman et al. (2018), mainly with respect
### the trade-off between accuracy and power analysis, since a p-value correction for 
### multiple tests (Sidak correction) is performed.

candidates <- listw.candidates(coord = ll, nb = c("del", "gab", "rel", "mst",
  			"pcnm", "dnear"), weights = c("binary", "flin", "fup", "fdown"))
names(candidates)                              
(nbw<-length(candidates)) ### Number of spatial weighting matrices generated

### Optimization the selection of a subset of SWM among the candidates generated above,
### using the corrected significance threshold calculated ("MIR"):
mod <- lm(Y ~ treat, data=matrix)
summary(mod)
res <- residuals(mod)
(W_sel_MIR <- listw.select(res, candidates, MEM.autocor = "positive", method = "MIR",
                           p.adjust = TRUE, MEM.all = FALSE, nperm = 999)) 
save.image()

### Some characteristics of the best spatial model:
# Best SWM:
W_sel_MIR$best.id
(names.sel<-names(W_sel_MIR$best.id))

# Retained object for further analysis (optional)
#SWM.selected<-W_sel_MIR$best$MEM.select
#class(SWM.selected)
#dim(SWM.selected)
#write.table(SWM.selected,"SWMselected.csv")

# Write the selected MEMs to a new object
#spatial.red <- as.matrix(SWM.selected)
#class(spatial.red)

### Testing the significance of the treatment, discouting the effect of 
### selected MEMs:
matrix2 <- cbind(matrix, W_sel_MIR$best$MEM.select)
head(matrix2)
mod_final <- lm(Y ~ treat + ., data = as.data.frame(matrix2))
anova<-aov(mod_final)
summary(anova)

### In case of significance, you should perform post hoc tests (Tukey etc.).
TukeyHSD(anova, "treat", ordered = TRUE)
plot(TukeyHSD(anova, "treat"))

# Important: don't forget to verify also the residuals normality and homogeneity of variances assumptions, since ANOVA is a parametric test!

save.image()
# End.

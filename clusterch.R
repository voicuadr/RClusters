########################################
#Downloading R packages
source("https://bioconductor.org/biocLite.R")
biocLite("ChemmineR")
library(rcdk)#
library(chemometrics)#
library(rJava)
library(ChemmineR) 
library(cluster)
library(rgl)
library(ggplot2)
library(vegan)
library(factoextra)
library(fingerprint)
library(fmcsR)
source("https://bioconductor.org/biocLite.R")
biocLite("fmcsR")
library(NbClust )
library(iqspr)
library(ggplot2)
library(gridExtra)
library(fpc)
###############################################
#Reading and visualising the methotrexate molecule in SMILES format 
mol<- parse.smiles('CCOC(=O)C1=C(CN=C1C)\\N=N\\C1=C(O)C=CC2=C1C=C(O)C=C2',kekulise=TRUE)[[1]]
mol
view.molecule.2d(mol)

#Reading and visualising the methotrexate (CMP1) in sdf format 
CMP1 <- load.molecules( c('CMP1.sdf') )
 view.molecule.2d(CMP1[[1]])

#########################################################
#General properties of methotrexate molecule
cat('No. of atoms =', length(atoms), '\n')
cat('No. of bonds =', length(bonds), '\n')
atoms <- get.atoms((CMP1[[1]]))
atoms
bonds <- get.bonds((CMP1[[1]]))
bonds
coords <- get.point2d(atoms[[1]])
coords
coords <- do.call('rbind'  , lapply(atoms, get.point2d))
coords
############################################################
#Descriptors categories in rcdk package:
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
descNames
dc <- get.desc.categories()
dc
descriptors = get.desc.names(type="all")
descriptors
#"Constitutional" descriptors 
dn <- get.desc.names(dc[2])
dn

#Calculus of a descriptor - 14 "AlogP"
aDesc <- eval.desc(CMP1, dn[14])
aDesc

# Topological Polar Surface Area, xlogp, alogp, total charge
get.tpsa(mol) 
get.xlogp(mol)
get.alogp(mol)
get.total.charge(mol)

#Calculus of all the descriptors  -1
allDescs <- eval.desc(CMP1, dn)
allDescs
require(rcdk)
#Calculus of all the descriptors  -2
drug.mols <- load.molecules(molfiles="CMP1.sdf")
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
drug.descs <- eval.desc(drug.mols, descNames, verbose=T)
drug.descs
descs <- eval.desc(mols, descNames)

##########################################################
#Calculs of the molecular fingerprint for methotrexate (CMP1)-maccs

fps <- get.fingerprint(CMP1[[1]], type='maccs')
fps
#Calculs of the molecular fingerprint for methotrexate (CMP1)-extended

fps <- get.fingerprint(CMP1[[1]], type='extended')
fps
## Reading and visualizing the set of 24 molecules
mols <- load.molecules(c('CMP1.sdf', 'CMP2.sdf', 'CMP3.sdf','CMP4.sdf', 'CMP5.sdf', 'CMP6.sdf', 'CMP7.sdf', 'CMP8.sdf', 'CMP9.sdf', 'CMP10.sdf', 'CMP11.sdf', 'CMP12.sdf', 'CMP13.sdf', 'CMP14.sdf', 'CMP15.sdf','CMP16.sdf', 'CMP17.sdf','CMP18.sdf', 'CMP19.sdf','CMP20.sdf','CMP21.sdf','CMP22.sdf', 'CMP23.sdf','CMP24.sdf'))

view.molecule.2d(mols, ncol = 4,width = 200, height = 200, depictor = NULL, type="isomeric")

view.molecule.2d(mols)
mols


#The molecular fingerprints for the set of molecules 
fps <- lapply(mols, get.fingerprint, type='extended')
fps
fp.sim <- fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim
fp.dist



#Identification of molecules located at a certain distance from the target
query.fp<-get.fingerprint(CMP1[[1]], type = 'maccs')
target.mols <-mols
target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
target.fps
sims <- data.frame(sim=do.call(rbind, lapply(target.fps,
                                             fingerprint::distance,
                                             fp2=query.fp, method='tanimoto')))
subset(sims, sim >= 0.5)
hits <- which(sims >= 0.5)
hits


#Distances between CMP1 and the rest of molecules
query.mol<-load.molecules( c('CMP1.sdf') )
query.mol
target.mols<-mols
target.mols
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim
fp.dist


####################################################################

#Assessment of the optimal number of clusters using NbClust R package
fviz_nbclust(fp.dist, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)
##############################################

#Hierarchical clustering with hclust using Ward's method

d <- dist(fp.dist, method = "euclidean")
res.hc <- hclust(d, method = "ward.D2" )
grp <- cutree(res.hc, k = 3)
plot(res.hc, cex = 0.6) 
rect.hclust(res.hc, k = 3, border = 2:5) 

#K-means clustering
fviz_nbclust(fp.dist, method = "gap_stat")
km.res <- kmeans(fp.dist, 3, nstart = 5)
km.res
fviz_cluster(km.res, data = fp.dist, ellipse.type = "convex")
theme_minimal()

####################################
#Clustering analysis
#Clusters statistics for K-mean

silinfo <- km.res$silinfo
names(silinfo)
km_stats <- cluster.stats(fp.dist,  km.res$cluster)
km_stats

# Silhouette coefficient of observations
library("cluster")
sil <- silhouette(km.res$cluster, dist(fp.dist))
head(sil[, 1:3], 10)
plot(sil, main ="Silhouette plot - K-means")

#####################################

#Visualisation of the molecules from cluster

molsc <- load.molecules(c('CMP1.sdf', 'CMP2.sdf', 'CMP3.sdf', 'CMP5.sdf', 'CMP6.sdf', 'CMP7.sdf', 'CMP8.sdf'))
view.molecule.2d(molsc)

##############################


















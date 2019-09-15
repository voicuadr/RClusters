########################################
#Incarcarea pachetelor necesare
library(rcdk)#
library(chemometrics)#
library(rJava)
library(ChemmineR) 
library(cluster)
library(rgl)
library(ggplot2)
library(vegan)
library(factoextra)
###############################################
#Citirea si vizualizarea moleculei"drug candidate"in format SMILES
mol<- parse.smiles('CCOC(=O)C1=C(CN=C1C)\\N=N\\C1=C(O)C=CC2=C1C=C(O)C=C2',kekulise=TRUE)[[1]]
#convert.implicit.to.explicit(mol)
mol
view.molecule.2d(mol)

#Citirea si vizualizarea moleculei de methotrexate in format sdf
meth <- load.molecules( c('meth.sdf') )
view.molecule.2d(meth[[1]])
#Vizualizarea moleculei informat SMILES
get.smiles(meth[[1]])

#########################################################
#Proprietatile generale ale moleculei
cat('No. of atoms =', length(atoms), '\n')
cat('No. of bonds =', length(bonds), '\n')
atoms <- get.atoms((meth[[1]]))
atoms
bonds <- get.bonds((meth[[1]]))
bonds
coords <- get.point2d(atoms[[1]])
coords
coords <- do.call('rbind'  , lapply(atoms, get.point2d))
coords
############################################################
#Categoriile de descriptori din pachetul rcdk:
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
descNames
dc <- get.desc.categories()
dc
descriptors = get.desc.names(type="all")
descriptors
#Descriptorii din categoria 2  "constitutional"
dn <- get.desc.names(dc[2])
dn

#Calculul unui descriptor - 14"AlogP"
aDesc <- eval.desc(meth, dn[14])
aDesc
#Calcularea tuturor descriptorilor
allDescs <- eval.desc(meth, dn)
allDescs
require(rcdk)

drug.mols <- load.molecules(molfiles="meth.sdf")
descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
drug.descs <- eval.desc(drug.mols, descNames, verbose=T)
drug.descs
descs <- eval.desc(mols, descNames)

get.tpsa(mol) # Topological Polar Surface Area
get.xlogp(mol)
get.alogp(mol)
get.total.charge(mol)
##########################################################
#Calculul amprentei moleculare pentru molecula candidat

fps <- get.fingerprint(meth[[1]], type='maccs')
fps
#Pentru compusi aromatici se prefera "extended"
fps <- get.fingerprint(meth[[1]], type='extended')
fps
## Citirea si vizualizarea unui set de 24 molecule
mols <- load.molecules(c('meth.sdf', '41.sdf', '42.sdf','43.sdf', '45.sdf', '57.sdf', 
                         '58.sdf', '60.sdf', '61.sdf', 'cfr2.sdf', 'CF20.sdf', 'CF22.sdf', 'CF23.sdf', 'CF24.sdf', 'CF29.sdf','CF30.sdf', 'CF31.sdf','CP2.sdf', 'CP3.sdf','CP5.sdf','CP6.sdf','CP7.sdf', 'CP8.sdf','CP9.sdf'))
view.molecule.2d(mols)
mols
#Obtinerea amprentelor moleculare pentru setul de molecule
fps <- lapply(mols, get.fingerprint, type='extended')
fps
fp.sim <- fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim
fp.dist

#Identificarea moleculelor situate la o anumita distanta fata de target 
#target.mols <-mols
#query.fp <- get.fingerprint(query.mols, type='maccs')
#target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
#sims <- unlist(lapply(target.fps, distance, fp2=query.fp, method='tanimoto'))
sims
hits <- which(sims > 0.4)
hits
query.mol<-load.molecules( c('meth.sdf') )
query.mol
target.mols<-mols
target.mols
#query.mol<- parse.smiles('CCOC(=O)C1=C(CN=C1C)\\N=N\\C1=C(O)C=CC2=C1C=C(O)C=C2')[[1]]
#smiles <- c('CCOC(=O)C1=C(CN=C1C)\\N=N\\C1=C(O)C=CC2=C1C=C(O)C=C2','CCOC(=O)C1=C2N=NC3=C(C=CC4=C3C=C(O)C=C4)N2N=C1C', 'CCOC(=O)C1=NNC(\\N=N\\C2=C(O)C=CC3=C2C2=C(C=C3)N3N=C(C)C(C(=O)CC)=C3N=N2)=C1C', 'CCC(=O)C1=C2N=NC3=C(OC)C=C(OC)C=C3N2N=C1C' , 'CCCOC(=O)C1=C(O)N2N=C(C)C(C(=O)OCC)=C2N=N1', 'CCOC(=O)C1=C2N=NC(C#N)=C(O)N2N=C1C', 'CCCCCCOCC1=C(O)N2N=C(C)C(C(=O)OCC)=C2N=N1', 'CCOC(=O)C1=C2N=NC(Cl)C(=O)N2N=C1C')
#target.mols <- parse.smiles(smiles)
#query.fp <- get.fingerprint(query.mol[[1]], type='maccs')
#target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
#sims <- unlist(lapply(target.fps, distance, fp2=query.fp, method='tanimoto'))
#sims
#hits <- which(sims > 0.5)
#hits

####################################################################

#Determinarea nr optim de clustere
fviz_nbclust(fp.dist, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

#Cluster dendograma 1
res <- hcut(fp.dist, k = 5, stand = TRUE)
res
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))
rect.hclust(H.fit, k=3, border="red")



#Dendrograma 3

d <- dist(fp.dist, method = "euclidean")
# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )

# Cut tree into 4 groups
grp <- cutree(res.hc, k = 3)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 3, border = 2:5) # add rectangle
###


#Cluster 2
wssplot(fp.dist)
km <- kmeans(fp.dist,3,1000) # Kmeans Cluster with 3 centers and iterations =1000
km
cmd<-cmdscale(fp.dist)
cmd
cols <- c("steelblue", "darkred", "darkgreen", "pink","green")
plot(cmd,type="n")
groups <- levels(factor(km$cluster))
for(i in seq_along(groups))
{
  points(cmd[factor(km$cluster) == groups[i], ], col = cols[i], pch = 16)
}
ordispider(cmd, factor(km$cluster), label = TRUE)
ordihull(cmd, factor(km$cluster), lty = "dotted")

#cluster 4

#Embedding using principal components
pc<-prcomp(fp.dist)
pc
plot(pc$x[,1], pc$x[,2],col=km$cluster,pch=16)
pc<-cbind(pc$x[,1], pc$x[,2])
ordispider(pc, factor(km$cluster), label = TRUE)
ordihull(fp.dist, factor(km$cluster), lty = "dotted")

######################
#Cluster poligonal 3

fviz_nbclust(fp.dist, method = "gap_stat")
km.res <- kmeans(fp.dist, 3, nstart = 5)
km.res
# Visualize
fviz_cluster(km.res, data = fp.dist, ellipse.type = "convex")
theme_minimal()

#Cluster
library(cluster)
fit <- kmeans(fp.dist, 5)
fit
clusplot(fp.dist, fit$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0)



#Cluster 3D
cmd3d<-cmdscale(fp.dist,k=3)
pc<-prcomp(fp.dist)
pc
pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
plot3d(pc3d, col = km$cluster,type="s",size=1,scale=0.2)
plot3d(cmd3d, col = km$cluster,type="s",size=1,scale=0.2)


#########












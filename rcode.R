#Specifying the package to work with
# the two files 2gaq.dcd as well as 2gaq.pdb are in the same directory as the R code
# Before running this code make to run getwd() and if it didn't specify the directory that
# your file is in run the command setwd('path to folder')
#Some commands are run between the code, but were not included
#The following code is only needed to duplicate the results that we obtained

require(bio3d)
dcdfile <- "2gaq.dcd"
pdbfile <- "2gaq.pdb"
dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

#Select all C-alpha atoms for trajectory frame superposition
#ca.inds stands for Carbon-alpha indices (I did not know that before hand)
ca.inds <- atom.select(pdb, elety="CA")
#Superposing all frames
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
	           fixed.inds=ca.inds$xyz,
	           mobile.inds=ca.inds$xyz)

##Root mean Square Deviation (RMSD)
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram",
     			xlab="RMSD")
lines(density(rd), col="gray", lwd=3)

##Root mean Square Fluctuations (RMSF)
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")

##Principal Component Analysis
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)))
#Perform Clustering
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)
#Examine the contribution of each residue to pc1
plot.bio3d(pc$au[,1], ylab="PC1 (A)", 
				xlab="Residue Position",
				typ="l")
points(pc$au[,2], typ="l",col="blue")
#To visualize the PCAs
p1 <- mktrj.pca(pc,pc=1,b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2, b=pc$au[,2], file="pc2.pdb")
p3 <- mktrj.pca(pc, pc=3, b=pc$au[,3], file="pc3.pdb")
p4 <- mktrj.pca(pc, pc=4, b=pc$au[,4], file="pc4.pdb")
p5 <- mktrj.pca(pc, pc=5, b=pc$au[,5], file="pc5.pdb")

##Cross Correlation Analysis
cij <- dccm(xyz[,ca.inds$xyz])
plot(cij)
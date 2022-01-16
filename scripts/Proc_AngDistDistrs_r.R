subj=commandArgs(trailingOnly=TRUE)
print(subj)

library(diptest)

# read in ang distrs for this subj
angdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))
angdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_R.csv'))
# and for group
gangdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked_L.csv'))
gangdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked_R.csv'))

# process full-brain-merged distributions
meltedL<-unlist(data.frame(angdL))
meltedR<-unlist(data.frame(angdR))
melted=c(meltedL,meltedR)
# and for group
gmeltedL<-unlist(data.frame(gangdL))
gmeltedR<-unlist(data.frame(gangdR))
gmelted=c(gmeltedL,gmeltedR)

# now we are working across all faces across both hemispheres

# divisor for binning 360 total degrees
divs=180/24

# initialize shape matrix
shapeMatrix<-array(0,dim=c(24,2))
gshapeMatrix<-array(0,dim=c(24,2))
# saveout x,y shape (24 bins, full distr)
shapeMatrix[,1]<-seq(1,24)
gshapeMatrix[,1]<-seq(1,24)

# discr. to 24 bins to eval. which bin is best represented
for (b in 1:24){
	range=c(((b-1)*divs),(b*divs))
	# set y coord to density of this bin
	shapeMatrix[b,2]=length(melted[melted < range[2] & melted > range[1]])
	gshapeMatrix[b,2]=length(gmelted[gmelted < range[2] & gmelted > range[1]])
}
 
# dip test statistic
DipT=dip.test(melted)
gDipT=dip.test(gmelted)
print(DipT)
print(gDipT)

# save shape matrix (subjs by 24 by 2)
saveRDS(shapeMatrix,paste0('~/results/PWs/Proced/',subj,'/',subj,'_ShapeMatrix.rds'))
saveRDS(gshapeMatrix,paste0('~/results/PWs/Proced/',subj,'/',subj,'_gShapeMatrix.rds'))
saveRDS(DipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_DipTest.rds'))
saveRDS(gDipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_gDipTest.rds'))

### dip test each face
# initialize vecs
dipvecL=array(0,length(angdL))
dipvecR=array(0,length(angdR))
gdipvecL=array(0,length(gangdL))
gdipvecR=array(0,length(gangdR))

# loop over left
for (f in 1:length(angdL)){
	dipSt=dip.test(angdL[,f])
	dipvecL[f]=dipSt$statistic
}
for (f in 1:length(gangdL)){
	gdipSt=dip.test(gangdL[,f])
	gdipvecL[f]=gdipSt$statistic
}
# loop over right
for (f in 1:length(angdR)){
        dipSt=dip.test(angdR[,f])
        gdipSt=dip.test(gangdR[,f])
}
for (f in 1:length(gangdR)){
        dipvecR[f]=dipSt$statistic
        gdipvecR[f]=gdipSt$statistic
}
# saveout
saveRDS(dipvecL,paste0('~/results/PWs/Proced/',subj,'/',subj,'_LVerts_DipTest.rds'))
saveRDS(gdipvecL,paste0('~/results/PWs/Proced/',subj,'/',subj,'_LVerts_gDipTest.rds'))
saveRDS(dipvecR,paste0('~/results/PWs/Proced/',subj,'/',subj,'_RVerts_DipTest.rds'))
saveRDS(gdipvecR,paste0('~/results/PWs/Proced/',subj,'/',subj,'_RVerts_gDipTest.rds'))


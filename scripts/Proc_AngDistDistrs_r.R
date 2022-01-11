subj=commandArgs(trailingOnly=TRUE)
print(subj)

library(diptest)

# read in ang distrs for this subj
angdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))
angdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_R.csv'))

# process full-brain-merged distributions
meltedL<-unlist(data.frame(angdL))
meltedR<-unlist(data.frame(angdR))
melted=c(meltedL,meltedR)
absmelted=abs(melted)
# now we are working across all faces across both hemispheres

# divisor for binning 360 total degrees
divs=360/24
divsAbs=180/24

# initialize shape matrix
shapeMatrix<-array(0,dim=c(24,2))
shapeMatrixAbs<-array(0,dim=c(24,2))
# saveout x,y shape (24 bins, full distr)
shapeMatrix[,1]<-seq(1,24)
shapeMatrixAbs[,1]<-seq(1,24)
# discr. to 24 bins to eval. which bin is best represented
for (b in 1:24){
	range=c((((b-1)*divs)-180),((b*divs)-180))
	rangeAbs=c((((b-1)*divsAbs)),((b*divsAbs)))
	# set y coord to density of this bin
	shapeMatrix[b,2]=length(melted[melted < range[2] & melted > range[1]])
	shapeMatrixAbs[b,2]=length(absmelted[absmelted < rangeAbs[2] & absmelted > rangeAbs[1]])
}
 
# dip test statistic
DipT=dip.test(absmelted)
print(DipT)

# save shape matrix (subjs by 24 by 2)
saveRDS(shapeMatrix,paste0('~/results/PWs/Proced/',subj,'/',subj,'_ShapeMatrix.rds'))
saveRDS(shapeMatrixAbs,paste0('~/results/PWs/Proced/',subj,'/',subj,'_AbsShapeMatrix.rds'))
saveRDS(DipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_DipTest.rds'))

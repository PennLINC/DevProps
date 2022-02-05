subj=commandArgs(trailingOnly=TRUE)
print(subj)

library(diptest)

# read in ang distrs for this subj
gangdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_c_L.csv'))
gangdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_c_R.csv'))

# process full-brain-merged distributions
gmeltedL<-unlist(data.frame(gangdL))
gmeltedR<-unlist(data.frame(gangdR))
gmelted=c(gmeltedL,gmeltedR)

# now we are working across all faces across both hemispheres

# divisor for binning 360 total degrees
divs=180/24

# initialize shape matrix
gshapeMatrix<-array(0,dim=c(24,2))
# saveout x,y shape (24 bins, full distr)
gshapeMatrix[,1]<-seq(1,24)

# discr. to 24 bins to eval. which bin is best represented
for (b in 1:24){
	range=c(((b-1)*divs),(b*divs))
	# set y coord to density of this bin
	gshapeMatrix[b,2]=length(gmelted[gmelted < range[2] & gmelted > range[1]])
}
 
# dip test statistic
gDipT=dip.test(gmelted)
print(gDipT)

# save shape matrix (subjs by 24 by 2)
saveRDS(gshapeMatrix,paste0('~/results/PWs/Proced/',subj,'/',subj,'_gShapeMatrix_c.rds'))
saveRDS(gDipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_gDipTest_c.rds'))

### dip test each face
# initialize vecs
gdipvecL=array(0,length(gangdL))
gdipvecR=array(0,length(gangdR))

# loop over left
for (f in 1:length(gangdL)){
	gdipSt=dip.test(gangdL[,f])
	gdipvecL[f]=gdipSt$statistic
}

# loop over right
for (f in 1:length(gangdR)){
        gdipSt=dip.test(gangdR[,f])
        gdipvecR[f]=gdipSt$statistic
}
# saveout
saveRDS(gdipvecL,paste0('~/results/PWs/Proced/',subj,'/',subj,'_LVerts_gDipTest_c.rds'))
saveRDS(gdipvecR,paste0('~/results/PWs/Proced/',subj,'/',subj,'_RVerts_gDipTest_c.rds'))


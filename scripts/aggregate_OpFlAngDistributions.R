# aggregate angular distributions, run shape pca on them

# read subject list
subjList=read.delim('~/PWs/hcpd_subj_list.txt')

# initialize shape matrix
shapeMatrix<-array(0,dim=c(24,2,dim(subjList)[1]))

# divisor for 24/360 possible degrees
divs=360/24

for (s in 1:dim(subjList)[1]){
  subj=subjList[s,1]
  print(s)
  # if this subject ran, add it to the shape matrix 
  if(file.exists(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))){
 	 # load in data
 	 angdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))
 	 angdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_R.csv'))
 	 # melt
 	 meltedL<-unlist(data.frame(angdL))
 	 meltedR<-unlist(data.frame(angdR))
 	 melted=c(meltedL,meltedR)
 	 # extract shape
 	 # 1-24 as x-axis (24 bins representing angular distance from bottom up prop)
 	 shapeMatrix[,1,s]<-seq(1,24)
 	 # discr. to 24 bins to eval. which bin is best represented
 	 for (b in 1:24){
 	   range=c((((b-1)*divs)-180),((b*divs)-180))
 	   # set y coord to density of this bin
 	   shapeMatrix[b,2,s]=length(melted[melted < range[2] & melted > range[1]])
 	 }
 	 # normalize size of [,2,] to 1:24 range
 	 countRange=range(shapeMatrix[,2,s])
 	 countRangeDif=countRange[2]-countRange[1]
 	 countRangeDivisor=countRangeDif/24
 	 shapeMatrix[,2,s]=shapeMatrix[,2,s]/countRangeDivisor
}
}

# save shape matrix (subjs by 24 by 2)
saveRDS(shapeMatrix,'~/results/PWs/shapeMat.rds')



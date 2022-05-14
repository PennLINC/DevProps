.libPaths('~/R_mgcv')
library(pracma)

# load in yeo7 memberships
y_L=read.csv('~/data/yeo7FaceBooleans_L.csv')
y_R=read.csv('~/data/yeo7FaceBooleans_R.csv')

# try looping over every subject and getting every angular distance from every face
# load in subj list
subjList=read.table('~/PWs/rs_subs.csv')

# for each yeo7 networks
for (y in 1:7){
	# initialize vector of all angular distances
	valuesVec=rep(0,19)
	# load angDist iteratively
	for (s in 1:dim(subjList)[1]){
	    print(s)
	    subj=subjList[s,1]
   	    ResFp=paste0('~/results/PWs/Proced/',subj,'/')
	    # if output exists
	    if (file.exists(paste0(ResFp,subj,'_gAngDist_Masked4_L.csv'))) {
	      AngDists_L<-read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_L.csv'))
	      AngDistsmat=t(AngDists_L[y_L[,y]])
	      histoVers<-histc(AngDistsmat,edges=seq(0,180,10))
	      valuesVec=valuesVec+rowSums(histoVers$cnt)
	      # and right hemi
	      AngDists_R<-read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_R.csv'))
	      AngDistsmat=t(AngDists_R[y_R[,y]])
	      histoVers<-histc(AngDistsmat,edges=seq(0,180,10))
	      valuesVec=valuesVec+rowSums(histoVers$cnt)
	    }
	}
	# convert to dataframe
	#fn
	fn=paste0('~/data/y_',y,'_allFacesAllSubjs.csv')
	# write csv
	write.csv(valuesVec,fn)
# end for yeo7 loop
}


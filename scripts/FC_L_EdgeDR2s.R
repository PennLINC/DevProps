iter=commandArgs(trailingOnly=TRUE)
.libPaths('~/R_mgcv')
library(mgcv)
library(ppcor)

# functions for delta r ^ 2
# difference in R2
DeltaR2EstVec<-function(x){
  
  # relevant df
  scaledf<-data.frame(cbind(as.numeric(masterdf$interview_age),as.numeric(masterdf$sex),masterdf$FD,x))
  colnames(scaledf)<-c('Age','Sex','Motion','varofint')
  
  # no-age model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Sex+Motion,data=scaledf)
  noAgeSum<-summary(noAgeGam)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=scaledf)
  AgeSum<-summary(AgeGam)
  
  dif<-AgeSum$r.sq-noAgeSum$r.sq
  
  # partial spearmans to extract age relation (for direction)
  pspear=pcor(scaledf,method='spearman')$estimate
  corest<-pspear[4]
  if(corest<0){
    dif=dif*-1
  }
  
  return(dif)
  
}

# same thing but returning chisq test sig. output for FDR correction instead of hard difference
DeltaPEstVec<-function(x){
  
  # relevant df
  scaledf<-data.frame(cbind(as.numeric(masterdf$interview_age),as.numeric(masterdf$sex),masterdf$FD,x))
  colnames(scaledf)<-c('Age','Sex','Motion','varofint')
  
  # no-age model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Sex+Motion,data=scaledf)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=scaledf)
  
  # test of dif with anova.gam
  anovaRes<-anova.gam(noAgeGam,AgeGam,test='Chisq')
  anovaP<-anovaRes$`Pr(>Chi)`
  anovaP2<-unlist(anovaP)
  return(anovaP2[2])
  
}

### load in subj info
# load in subj list
subjList=read.delim('/cbica/projects/pinesParcels/PWs/hcpd_subj_list.txt',header=F)

# load in ages
demo=read.csv('/cbica/projects/pinesParcels/PWs/hcpd_demographics.csv')
# convert to used naming convention
demo$SubjID<-gsub('HCD','sub-',demo$src_subject_id)

# load in FD
FD_TRs=read.csv('/cbica/projects/pinesParcels/PWs/Subj_FD_RemTRs.csv')
colnames(FD_TRs)[1]<-'SubjID'
colnames(FD_TRs)[2]<-'FD'
colnames(FD_TRs)[3]<-'RemainingTRs'

# merge by subjID
mergeddf<-merge(demo,FD_TRs,by='SubjID')

# exclude subjects with less than 600 TRs remaining (in rest)
inclusionVec<-mergeddf$RemainingTRs>600
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
mergeddf<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(mergeddf)[1]
print(remainingSubjs)

# make an edge matrix (subjs x number edges), use one participants edge file as template
subj=mergeddf$SubjID[1]
templatefileFP=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt', as.character(iter), '.rds')
templatefile=readRDS(templatefileFP)
edgedf=matrix(0,remainingSubjs,length(templatefile))
edgedf=data.frame(edgedf)
edgedf$SubjID=mergeddf$SubjID

### load in FC matrices (vectorized) for each participant
for (s in 1:remainingSubjs){
	subj=edgedf$SubjID[s]
	# this iterations chunk of edges
	fp=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt', as.character(iter), '.rds')
	subjEdges=readRDS(fp)
	# little extra text in indexing to spare last column (subj ID)
	edgedf[s,1:(ncol(edgedf)-1)]=as.numeric(subjEdges)	
	print(s)
}

# merge dfs
masterdf<-merge(mergeddf,edgedf,by='SubjID')

# initialize output df (dr2, p)
outdf=data.frame(matrix(0,2,length(templatefile)))
### loop over each edge
for (e in 1:length(templatefile)){
	# extract DR2 (starts at column after OG merged df)
	outdf[1,e]=DeltaR2EstVec(masterdf[,(ncol(mergeddf)+e)])	
	# extract Pvalue
	outdf[2,e]=DeltaPEstVec(masterdf[,(ncol(mergeddf)+e)])
	print(e)
	print(iter)
}

### print out dr2 vector
DR2VecName=paste('/cbica/projects/pinesParcels/results/PWs/Age_EdgeLeve_FC_L_DR2_',iter,'.rds',sep='')
saveRDS(outdf,DR2VecName)


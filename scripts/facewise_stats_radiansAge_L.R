# so gr8ful i can do this w/o pmacs: see https://pennlinc.github.io/docs/cubic#using-rr-studio-and-installation-of-r-packages to customize r packages as needed
.libPaths('~/R_mgcv')
# mgcv needed for gams
library(mgcv)
# ppcor for directionality detection
library(ppcor)

# difference in R2
DeltaR2EstVec<-function(x){
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$interview_age),as.numeric(df$sex),df$FD,x))
  colnames(testdf)<-c('Age','Sex','Motion','varofint')
  # no-age model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Sex+Motion,data=testdf)
  noAgeSum<-summary(noAgeGam)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=testdf)
  AgeSum<-summary(AgeGam)
  dif<-AgeSum$r.sq-noAgeSum$r.sq
  # partial spearmans to extract age relation (for direction)
  pspear=pcor(testdf,method='spearman')$estimate
  corest<-pspear[4]
  if(corest<0){
    dif=dif*-1
  }
  return(dif) 
}

# Next, to derive the statistical significance of observed age effects, we need to test if the two models (one with an age term, one without) are significantly different. We use an ANOVA for this procedure. These p-values will eventually be FDR-corrected.

# chisq test sig. output
DeltaPEstVec<-function(x){
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$interview_age),as.numeric(df$sex),df$FD,x))  
  colnames(testdf)<-c('Age','Sex','Motion','varofint')
  # no-age model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Sex+Motion,data=testdf)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=testdf)  
  # test of dif with anova.gam
  anovaRes<-anova.gam(noAgeGam,AgeGam,test='Chisq')
  anovaP<-anovaRes$`Pr(>Chi)`
  anovaP2<-unlist(anovaP)
  return(anovaP2[2])  
}

# difference in R2
DeltaR2EstVec_s<-function(x){
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$sex),as.numeric(df$interview_age),df$FD,x))
  colnames(testdf)<-c('Sex','Age','Motion','varofint')
  # no-SEX model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Motion+s(Age,k=4),data=testdf)
  noAgeSum<-summary(noAgeGam)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=testdf)
  AgeSum<-summary(AgeGam)
  dif<-AgeSum$r.sq-noAgeSum$r.sq
  # partial spearmans to extract age relation (for direction)
  pspear=pcor(testdf,method='spearman')$estimate
  corest<-pspear[4]
  if(corest<0){
    dif=dif*-1
  }
  return(dif)
}

# Next, to derive the statistical significance of observed age effects, we need to test if the two models (one with an age term, one without) are significantly different. We use an ANOVA for this procedure. These p-values will eventually be FDR-corrected.

# chisq test sig. output
DeltaPEstVec_s<-function(x){
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$sex),as.numeric(df$interview_age),df$FD,x))
  colnames(testdf)<-c('Sex','Age','Motion','varofint')
  # no-SEX model (segreg ~ sex + motion)
  noAgeGam<-gam(varofint~Motion+s(Age,k=4),data=testdf)
  # age-included model for measuring difference
  AgeGam<-gam(varofint~Sex+Motion+s(Age,k=4),data=testdf)
  # test of dif with anova.gam
  anovaRes<-anova.gam(noAgeGam,AgeGam,test='Chisq')
  anovaP<-anovaRes$`Pr(>Chi)`
  anovaP2<-unlist(anovaP)
  return(anovaP2[2])
}

### end functions - start actual script ###

# extract range of vertices to be covered in this run from VertBin
Lfaces=4589
Rfaces=4595
# load in subj list
subjList=read.delim('/cbica/projects/pinesParcels/PWs/hcpd_subj_list.txt')
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
# exclude subjects with less than 600 TRs remaining
inclusionVec<-mergeddf$RemainingTRs>600
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
df<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(df)[1]
print(remainingSubjs)
# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
# for plotting means
RadDist=rep(0,Lfaces)
# for plotting age effect sizes
RadDist_adr2=rep(0,Lfaces)
# sex
RadDist_sdr2=rep(0,Lfaces)
# for fdr-correcting age associations
RadDist_ap=rep(0,Lfaces)
RadDist_sp=rep(0,Lfaces)
# subjvec to run in parallel for even more confidence in merging
Subjvec=rep(0,remainingSubjs)
# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$RadDist=rep(0,remainingSubjs)
# for each face in this run's range
for (f in 1:Lfaces){
	print(f)
	# load in D and shapes iteratively
	for (s in 1:remainingSubjs){
	    subj=df$SubjID[s]
	    ResFP_L=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'AngDist90_L.csv')
	    # if output exists
	    if (file.exists(paste0(ResFP_L))) {
	      # load in dat data
	      Res=read.csv(ResFP_L)
		# extract TD resvec length
		df$RadDist[s]=Res[f,1]	
	    }
	}
	# extract mean
	RadDist[f]=mean(df$RadDist)
        # extract age dr2
        RadDist_adr2[f]=DeltaR2EstVec(df$RadDist)
        # extract age p
        RadDist_ap[f]=DeltaPEstVec(df$RadDist)
	# sex effects
	RadDist_sdr2[f]=DeltaR2EstVec_s(df$RadDist)
	RadDist_sp[f]=DeltaPEstVec_s(df$RadDist)
}
# saveout means
write.csv(RadDist,'~/results/PWs/MeanRadDist_L.csv',col.names=F,row.names=F,quote=F)
# saveout dr2s and ps - still needs to be merged with results from other hemi for MC correction
saveRDS(RadDist_adr2,paste0('/cbica/projects/pinesParcels/results/PWs/RadDist_adr2.rds'))
saveRDS(RadDist_sdr2,'/cbica/projects/pinesParcels/results/PWs/RadDist_sdr2.rds')
saveRDS(RadDist_ap,paste0('/cbica/projects/pinesParcels/results/PWs/RadDist_p.rds'))
saveRDS(RadDist_sp,'/cbica/projects/pinesParcels/results/PWs/RadDist_sp.rds')


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


# difference in R2 for EF
EFDeltaR2EstVec<-function(x){
  
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$F1),as.numeric(df$interview_age),as.numeric(as.factor(df$sex)),df$FD,x))
  colnames(testdf)<-c('EF','Age','Sex','Motion','varofint')
  
  # no-EF model (b.w. ~ sex + motion)
  noEFGam<-gam(varofint~Sex+Motion+s(Age,k=3),data=testdf)
  noEFSum<-summary(noEFGam)
  # EF-included model for measuring difference
  EFGam<-gam(varofint~EF+Sex+Motion+s(Age,k=3),data=testdf)
  EFSum<-summary(EFGam)
  
  dif<-EFSum$r.sq-noEFSum$r.sq
  
  # partial spearmans to extract EF relation (for direction)
  pspear=pcor(testdf,method='spearman')$estimate
  corest<-pspear[5]
  if(corest<0){
    dif=dif*-1
  }
  
  return(dif)
  
}

# same thing but returning chisq test sig. output for FDR correction instead of hard difference
EFDeltaPEstVec<-function(x){
  
  # relevant df
  testdf<-data.frame(cbind(as.numeric(df$F1),as.numeric(df$interview_age),as.numeric(as.factor(df$sex)),df$FD,x))
  colnames(testdf)<-c('EF','Age','Sex','Motion','varofint')
  
  # no-EF model (segreg ~ sex + motion)
  noEFGam<-gam(varofint~Sex+Motion+s(Age,k=3),data=testdf)
  # EF-included model for measuring difference
  EFGam<-gam(varofint~EF+Sex+Motion+s(Age,k=3),data=testdf)
  
  # test of dif with anova.gam
  anovaRes<-anova.gam(noEFGam,EFGam,test='Chisq')
  anovaP<-anovaRes$`Pr(>Chi)`
  anovaP2<-unlist(anovaP)
  return(anovaP2[2])
  
}

### end functions - start actual script ###

# extract range of vertices to be covered in this run from VertBin
Lfaces=4851
Rfaces=4842

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

# load in Cog
Cog=read.delim('/cbica/projects/pinesParcels/data/cogcomp01.txt')
Cog$SubjID=gsub('HCD','sub-',Cog$src_subject_id)
Cog$F1<-Cog$nih_crycogcomp_unadjusted

# merge by subjID
mergeddf<-merge(demo,FD_TRs,by='SubjID')
mergeddf<-merge(mergeddf,Cog,by='SubjID')
# exclude subjects with less than 600 TRs remaining
inclusionVec<-mergeddf$RemainingTRs>600
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
df<-mergeddf[inclusionVec,]
# a little extra exclusion bc some subjs have no cog Comp scores
df<-df[complete.cases(as.numeric(df$F1)),]
# some column renaming to account for redundancy
df$interview_age<-df$interview_age.x
df$sex<-df$sex.x

# get count of remaining subjects
remainingSubjs=dim(df)[1]
print(remainingSubjs)

# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
# for plotting means
TD_R=rep(0,Rfaces)
BU_R=rep(0,Rfaces)
BuProp=rep(0,Rfaces)
ThetasFromPG=rep(0,Rfaces)

# for plotting age effect sizes
TD_R_edr2=rep(0,Rfaces)
BU_R_edr2=rep(0,Rfaces)
BuProp_edr2=rep(0,Rfaces)
ThetasFromPG_edr2=rep(0,Rfaces)

# for fdr-correcting age associations
TD_R_ep=rep(0,Rfaces)
BU_R_ep=rep(0,Rfaces)
BuProp_ep=rep(0,Rfaces)
ThetasFromPG_ep=rep(0,Rfaces)

# subjvec to run in parallel for even more confidence in merging
Subjvec=rep(0,remainingSubjs)

# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$FaceBuProp=rep(0,remainingSubjs)
df$FaceBu_rv=rep(0,remainingSubjs)
df$FaceTd_rv=rep(0,remainingSubjs)
df$FaceThetaDist=rep(0,remainingSubjs)

# that leaves a column for face value, one for age, one for FD, one for sex. Remaining TRs controlled for via exclusion

# for each face in this run's range
for (f in 1:4842){
	print(f)
	# load in D and shapes iteratively
	for (s in 1:remainingSubjs){
	    subj=df$SubjID[s]
	    ResFP_L=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_BUTD_R.csv')
	    # if output exists
	    if (file.exists(paste0(ResFP_L))) {
	      # load in dat data
	      Res=read.csv(ResFP_L)
		# extract TD resvec length
		df$FaceBuProp[s]=Res[f,1]	
		# extract BU resvec length
		df$FaceBu_rv[s]=Res[f,2]
		# extract prop of BU TRs
		df$FaceTd_rv[s]=Res[f,3]
		# extract global resvec theta from gPGG
		df$FaceThetaDist[s]=Res[f,4]
	    }
	}
        # extract EF dr2
        TD_R_edr2[f]=EFDeltaR2EstVec(df$FaceTd_rv)
        # extract age dr2
        BU_R_edr2[f]=EFDeltaR2EstVec(df$FaceBu_rv)
        # extract age dr2
        BuProp_edr2[f]=EFDeltaR2EstVec(df$FaceBuProp)
        # you already know doe
        ThetasFromPG_edr2[f]=EFDeltaR2EstVec(df$FaceThetaDist)
	
	# extract EF p
        TD_R_ep[f]=EFDeltaPEstVec(df$FaceTd_rv)
        # extract age p
        BU_R_ep[f]=EFDeltaPEstVec(df$FaceBu_rv)
        # extract age p
        BuProp_ep[f]=DeltaPEstVec(df$FaceBuProp)
        # you already know doe
        ThetasFromPG_ep[f]=EFDeltaPEstVec(df$FaceThetaDist)

}

# saveout dr2s and ps - still needs to be merged with results from other hemi for MC correction
saveRDS(TD_R_edr2,paste0('/cbica/projects/pinesParcels/results/PWs/RTDL_cedr2.rds'))
saveRDS(BU_R_edr2,paste0('/cbica/projects/pinesParcels/results/PWs/RBUL_cedr2.rds'))
saveRDS(BuProp_edr2,paste0('/cbica/projects/pinesParcels/results/PWs/RBUProp_cedr2.rds'))
saveRDS(ThetasFromPG_edr2,paste0('/cbica/projects/pinesParcels/results/PWs/RThetasFromPG_cedr2.rds'))

saveRDS(TD_R_ep,paste0('/cbica/projects/pinesParcels/results/PWs/RTDL_cep.rds'))
saveRDS(BU_R_ep,paste0('/cbica/projects/pinesParcels/results/PWs/RBUL_cep.rds'))
saveRDS(BuProp_ep,paste0('/cbica/projects/pinesParcels/results/PWs/RBUProp_cep.rds'))
saveRDS(ThetasFromPG_ep,paste0('/cbica/projects/pinesParcels/results/PWs/RThetasFromPG_cep.rds'))


# so gr8ful i can do this w/o pmacs: see https://pennlinc.github.io/docs/cubic#using-rr-studio-and-installation-of-r-packages to customize r packages as needed
#.libPaths('~/R_mgcv')

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
# FD for carit runs
FD_TRs_c=read.csv('/cbica/projects/pinesParcels/PWs/Subj_FD_RemTRs_c.csv')
colnames(FD_TRs_c)[1]<-'SubjID'
colnames(FD_TRs_c)[2]<-'FDc'
colnames(FD_TRs_c)[3]<-'RemainingTRsc'

# merge by subjID
mergeddf<-merge(demo,FD_TRs,by='SubjID')
mergeddf<-merge(mergeddf,FD_TRs_c,by='SubjID')

# exclude subjects with less than 600 TRs remaining (in rest)
inclusionVec<-mergeddf$RemainingTRs>600
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
mergeddf<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(mergeddf)[1]
print(remainingSubjs)

# exclude subjects with less than 300 TRs remaining (in carit)
inclusionVec<-mergeddf$RemainingTRsc>300
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
df<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(df)[1]
print(remainingSubjs)

# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
T_L=rep(0,Lfaces)
T_R=rep(0,Rfaces)
p_L=rep(0,Lfaces)
p_R=rep(0,Rfaces)
dif_L=rep(0,Lfaces)
dif_R=rep(0,Rfaces)

# subjvec to run in parallel for even more confidence in merging
Subjvec=rep(0,remainingSubjs)

# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$BuProp=rep(0,remainingSubjs)
df$BuProp_c=rep(0,remainingSubjs)
# that leaves a column for face value, one for age, one for FD, one for sex. Remaining TRs controlled for via exclusion

# for each left hemi face
for (f in 1:Lfaces){
	print(f)
	# load in BU props iteratively
	for (s in 1:remainingSubjs){
	    subj=df$SubjID[s]
	    ResFP_L=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_BUTD_L.csv')
	    # load in dat data
	    Res=read.csv(ResFP_L)
    	    ResFP_L_c=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_BUTD_L_c.csv')
            # load in dat data
            Res_c=read.csv(ResFP_L_c)
	    # extract rest BuProp length
	    df$BuProp[s]=Res[f,1]	
	    # extract carit BuProp length
	    df$BuProp_c[s]=Res_c[f,1]
	}
	# t test em
	ttestres=t.test(df$BuProp,df$BuProp_c,paired=TRUE) 
	# extract t stat
	T_L[f]=ttestres$statistic
	p_L[f]=ttestres$p.value
	dif_L[f]=ttestres$estimate
}

# saveout dif
write.csv(dif_L,'~/results/PWs/r_vs_c_dif_L.csv',col.names=F,row.names=F,quote=F)

# saveout dr2s and ps - still needs to be merged with results from other hemi for MC correction
saveRDS(T_L,paste0('/cbica/projects/pinesParcels/results/PWs/r_vs_c_t_L.rds'))
saveRDS(p_L,paste0('/cbica/projects/pinesParcels/results/PWs/r_vs_c_p_L.rds'))


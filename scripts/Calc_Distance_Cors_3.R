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

# load in general euclidean distance and heirarchical distance matrices
HD_v_l=readRDS('~/data/HD_v_l.rds')
HD_v_r=readRDS('~/data/HD_v_r.rds')
ED_v_l=readRDS('~/results/PWs/FC/vert_EucD_left_vec.rds')
ED_v_r=readRDS('~/results/PWs/FC/vert_EucD_right_vec.rds')
#### Load in distances
HD_f_l=readRDS('~/data/HD_f_l.rds')
HD_f_r=readRDS('~/data/HD_f_r.rds')
#### Load in Euclidean distances
ED_f_l=readRDS('~/results/PWs/face_EucD_left_vec.rds')
ED_f_r=readRDS('~/results/PWs/face_EucD_right_vec.rds')

# load in medial wall mask for vertices
mw_l_verts=as.logical(read.csv('~/data/mw_boolean_l.csv'))
mw_r_verts=as.logical(read.csv('~/data/mw_boolean_r.csv'))

# initialize output vectors
mergeddf$EucxFc<-rep(99,388)
mergeddf$HDxFc<-rep(99,388)
mergeddf$EucxCFc<-rep(99,388)
mergeddf$HDxCFc<-rep(99,388)

# loop over every subj
for (s in 101:150){
	subj=mergeddf$SubjID[s]
	print(subj)
	# load in cfc
	cfcFP_l=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L.csv')
	cfcFP_r=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R.csv')
	cfc_l=read.csv(cfcFP_l)
	cfc_r=read.csv(cfcFP_r)
	# # upper tri
	cfc_l=cfc_l[upper.tri(cfc_l)]
	cfc_r=cfc_r[upper.tri(cfc_r)]
	# load in fc
	fcfp=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_FC.csv')
	fc=read.csv(fcfp,header=F)
	# convert to only vertices of interest
	# left
	fc_l=fc[1:10242,1:10242]
	# mask out medial wall
	fc_l=fc_l[mw_l_verts,mw_l_verts]
	fc_l=fc_l[upper.tri(fc_l)]
	# right
	fc_r=fc[10243:20484,10243:20484]
	# mask out medial wall
	fc_r=fc_r[mw_r_verts,mw_r_verts]
	fc_r=fc_r[upper.tri(fc_r)]
	# calc spearman's corrs for
	# euc x fc
	efc_s_l=cor.test(fc_l,ED_v_l,method='spearman')
	efc_s_r=cor.test(fc_r,ED_v_r,method='spearman')
	# hd x fc
	hdfc_s_l=cor.test(fc_l,HD_v_l,method='spearman')
        hdfc_s_r=cor.test(fc_r,HD_v_r,method='spearman')
	# euc x cfc
	ecfc_s_l=cor.test(cfc_l,ED_f_l,method='spearman')
        ecfc_s_r=cor.test(cfc_r,ED_f_r,method='spearman')
	# hd x cfc
	hdcfc_s_l=cor.test(cfc_l,HD_f_l,method='spearman')
        hdcfc_s_r=cor.test(cfc_r,HD_f_r,method='spearman')
	# plop into df, mean of both hemis
	mergeddf$EucxFc[s]=mean(c(efc_s_l$estimate,efc_s_r$estimate))
	mergeddf$HDxFc[s]=mean(c(hdfc_s_l$estimate,hdfc_s_r$estimate))
	mergeddf$EucxCFc[s]=mean(c(ecfc_s_l$estimate,ecfc_s_r$estimate))
        mergeddf$HDxCFc[s]=mean(c(hdcfc_s_l$estimate,hdcfc_s_r$estimate))	
}
# saveout dataframe
saveRDS(mergeddf[101:151,],'~/data/mergedDF_withDistanceCorrs_3.rds')

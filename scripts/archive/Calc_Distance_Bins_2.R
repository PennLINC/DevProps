.libPaths('~/R_mgcv')
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

# load in general euclidean distance distance matrices
ED_v_l=readRDS('~/results/PWs/FC/vert_EucD_left_vec.rds')
ED_v_r=readRDS('~/results/PWs/FC/vert_EucD_right_vec.rds')

#### Load in Euclidean distances
ED_f_l=readRDS('~/results/PWs/face_EucD_left_vec.rds')
ED_f_r=readRDS('~/results/PWs/face_EucD_right_vec.rds')

# ##### ENSURE IT NEEDS MASKING

# load in medial wall mask for vertices
mw_l_verts=as.logical(read.csv('~/data/mw_boolean_lfs4.csv'))
mw_r_verts=as.logical(read.csv('~/data/mw_boolean_rfs4.csv'))


# CREATE DISTANCE BINS (1:180 MM)
DistThreshs=seq(0,180,10)

##### OUTPUT WILL BE DISTANCE BIN BY SUBJ
outMat<-matrix(0,50,19)
outDf_fc<-data.frame(outMat)
outDf_cfc<-data.frame(outMat)

# initialize output vectors
#mergeddf$EucxFc<-rep(99,388)
#mergeddf$HDxFc<-rep(99,388)
#mergeddf$EucxCFc<-rep(99,388)
#mergeddf$HDxCFc<-rep(99,388)

# loop over every subj
for (s in 51:100){
	subj=mergeddf$SubjID[s]
	print(subj)
	age=mergeddf$interview_age[s]
	FD=mergeddf$FD[s]
	# load in cfc
	cfcFP_l=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L.csv')
	cfcFP_r=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R.csv')
	cfc_l=read.csv(cfcFP_l,header=F)
	cfc_r=read.csv(cfcFP_r,header=F)
	# # upper tri
	cfc_l=cfc_l[upper.tri(cfc_l)]
	cfc_r=cfc_r[upper.tri(cfc_r)]
	# load in fc
	fcfp=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_FCfs4.csv')
	fc=read.csv(fcfp,header=F)
	# convert to only vertices of interest
	# left
	fc_l=fc[1:2562,1:2562]
	fc_l=fc_l[mw_l_verts,mw_l_verts]
	fc_l=fc_l[upper.tri(fc_l)]
	# right
	fc_r=fc[2563:5124,2563:5124]
	# mask out medial wall
	fc_r=fc_r[mw_r_verts,mw_r_verts]
	fc_r=fc_r[upper.tri(fc_r)]
	# Make temporary dataframe to filter FC by distances
	tempDf_vl<-data.frame(fc_l,ED_v_l)
	tempDf_vr<-data.frame(fc_r,ED_v_r)
	tempDf_fl<-data.frame(cfc_l,ED_f_l)
	tempDf_fr<-data.frame(cfc_r,ED_f_r)
	# make colnames reasonable for filtering
	colnames(tempDf_vl)<-c('fc','ED')
	colnames(tempDf_vr)<-c('fc','ED')
	colnames(tempDf_fl)<-c('cfc','ED')
	colnames(tempDf_fr)<-c('cfc','ED')
	
	# FOR EACH DISTANCE BIN	
	for (db in 1:16){
		print(db)
		# Thresh Lower
		ThL=DistThreshs[db]
		# Thresh Higher
		ThH=DistThreshs[db+1]
		### FILTER OUT AVERAGE FC IN DISTANCE
		# Left
		FC_L=tempDf_vl[tempDf_vl$ED<ThH,]
		FC_L=FC_L[FC_L$ED>ThL,]		
		# Right
		FC_R=tempDf_vr[tempDf_vr$ED<ThH,]
                FC_R=FC_R[FC_R$ED>ThL,]
		### FILTER OUT AVERAGE CFC IN DISTANCE
		# Left
		CFC_L=tempDf_fl[tempDf_fl$ED<ThH,]
                CFC_L=CFC_L[CFC_L$ED>ThL,]
		# Right
		CFC_R=tempDf_fr[tempDf_fr$ED<ThH,]
                CFC_R=CFC_R[CFC_R$ED>ThL,]
		# average across hemis
		FC=mean(c(FC_L$fc,FC_R$fc))
		CFC=mean(c(CFC_L$cfc,CFC_R$cfc))
		# ENTER FC AND CFC INTO BIN VECTOR FOR THIS SUBJ
		outDf_fc[s,db]=FC
		outDf_cfc[s,db]=CFC
	}
	# END FOR EACH DISTANCE BIN	
	# enter subj ID
	outDf_fc[s,17]=subj
	outDf_cfc[s,17]=subj
	# enter age
	outDf_fc[s,18]=age
        outDf_cfc[s,18]=age
	# enter motion
	outDf_fc[s,19]=FD
        outDf_cfc[s,19]=FD
	# END FOR EACH SUBJ
}

# SAVEOUT BINNED FCS AND CFCS
saveRDS(outDf_fc,'~/data/outDf_fc_2.rds')
saveRDS(outDf_cfc,'~/data/outDf_cfc_2.rds')

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

# load in medial wall mask for vertices
mw_l_verts=as.logical(read.csv('~/data/mw_boolean_lfs4.csv'))
mw_r_verts=as.logical(read.csv('~/data/mw_boolean_rfs4.csv'))

# loop over every subj
for (s in 1:remainingSubjs){
        subj=mergeddf$SubjID[s]
        print(subj)
	# read cfc
        cfcFP_l=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L.csv')
        cfcFP_r=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R.csv')
        cfc_l=read.csv(cfcFP_l,header=F)
        cfc_r=read.csv(cfcFP_r,header=F)
        # # upper tri
        cfc_l=cfc_l[upper.tri(cfc_l)]
        cfc_r=cfc_r[upper.tri(cfc_r)]

	# save upper triangles out left
	saveRDS(cfc_l[1:1052718],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt1.rds'))
	saveRDS(cfc_l[1052719:2105435],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt2.rds'))
	saveRDS(cfc_l[2105436:3158152],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt3.rds'))
	saveRDS(cfc_l[3158153:4210869],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt4.rds'))
	saveRDS(cfc_l[4210870:5263586],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt5.rds'))
	saveRDS(cfc_l[5263587:6316303],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt6.rds'))
	saveRDS(cfc_l[6316304:7369020],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt7.rds'))
	saveRDS(cfc_l[7369021:8421737],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt8.rds'))
	saveRDS(cfc_l[8421738:9474454],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt9.rds'))
	saveRDS(cfc_l[9474455:length(cfc_l)],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_L_pt10.rds'))

	# save upper triangles out right
        saveRDS(cfc_r[1:1052718],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt1.rds'))
        saveRDS(cfc_r[1052719:2105435],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt2.rds'))
        saveRDS(cfc_r[2105436:3158152],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt3.rds'))
        saveRDS(cfc_r[3158153:4210869],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt4.rds'))
        saveRDS(cfc_r[4210870:5263586],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt5.rds'))
        saveRDS(cfc_r[5263587:6316303],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt6.rds'))
        saveRDS(cfc_r[6316304:7369020],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt7.rds'))
        saveRDS(cfc_r[7369021:8421737],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt8.rds'))
        saveRDS(cfc_r[8421738:9474454],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt9.rds'))
        saveRDS(cfc_r[9474455:length(cfc_r)],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_CircFC_R_pt10.rds'))

	# load in fc
        fcfp=paste0('/cbica/projects/pinesParcels/results/PWs/Proced/',subj,'/',subj,'_FCfs4.csv')
        fc=read.csv(fcfp,header=F)
        # convert to only vertices of interest
        # left
        fc_l=fc[1:2562,1:2562]
        # mask out medial wall
        fc_l=fc_l[mw_l_verts,mw_l_verts]
        fc_l=fc_l[upper.tri(fc_l)]
        # right
        fc_r=fc[2563:5124,2563:5124]
        # mask out medial wall
        fc_r=fc_r[mw_r_verts,mw_r_verts]
        fc_r=fc_r[upper.tri(fc_r)]

	# save out broken up version
	saveRDS(fc_l[1:542201],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt1.rds'))
        saveRDS(fc_l[542202:1084401],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt2.rds'))
        saveRDS(fc_l[1084402:1626601],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt3.rds'))
        saveRDS(fc_l[1626602:2168801],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt4.rds'))
        saveRDS(fc_l[2168802:length(fc_l)],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_L_pt5.rds'))

	# right
	saveRDS(fc_r[1:542201],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_R_pt1.rds'))
        saveRDS(fc_r[542202:1084401],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_R_pt2.rds'))
        saveRDS(fc_r[1084402:1626601],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_R_pt3.rds'))
        saveRDS(fc_r[1626602:2168801],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_R_pt4.rds'))
        saveRDS(fc_r[2168802:length(fc_l)],paste0('/cbica/projects/pinesParcels/results/PWs/Proced/', subj, '/', subj, '_FC_R_pt5.rds'))
}


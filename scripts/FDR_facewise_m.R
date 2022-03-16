# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold dr2 values accordingly

# load in facewise output from left hemi
LBP_m=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_mdr2.rds')
LBP_mp=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_mp.rds')

# and right hemi
RBP_m=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_mdr2.rds')
RBP_mp=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_mp.rds')

# combine each
BUs=c(LBP_m,RBP_m)
BUsP=c(LBP_mp,RBP_mp)

# fdr each
BUsP_f=p.adjust(BUsP,method='fdr')

# mask dr2s accordingly
BUsP_f[BUsP_f>0.05]=0

# uncombine: seperate vecs for sep. hemis
BU_L=BUsP[1:length(LBP_m)]
BU_R=BUsP[(length(LBP_m)+1):length(BUsP)]

# print out each for matlab friendly reading
write.table(BU_L,'~/results/PWs/FDRed_m_L.csv',col.names=F,row.names=F,quote=F)
write.table(BU_R,'~/results/PWs/FDRed_m_R.csv',col.names=F,row.names=F,quote=F)


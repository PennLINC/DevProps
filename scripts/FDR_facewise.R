# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold dr2 values accordingly

# load in facewise output from left hemi
LBU_L_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUL_adr2.rds')
LTD_L_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LTDL_adr2.rds')
LBuProp_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_adr2.rds')
LThetasFromPG_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LThetasFromPG_adr2.rds')

LBU_L_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUL_p.rds')
LTD_L_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LTDL_p.rds')
LBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_p.rds')
LThetasFromPG_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LThetasFromPG_p.rds')

# and right hemi
RBU_L_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUL_adr2.rds')
RTD_L_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RTDL_adr2.rds')
RBuProp_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_adr2.rds')
RThetasFromPG_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RThetasFromPG_adr2.rds')

RBU_L_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUL_p.rds')
RTD_L_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RTDL_p.rds')
RBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_p.rds')
RThetasFromPG_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RThetasFromPG_p.rds')

# combine each
BUdr2=c(LBU_L_adr2,RBU_L_adr2)
BUp=c(LBU_L_p,RBU_L_p)
TDdr2=c(LTD_L_adr2,RTD_L_adr2)
TDp=c(LTD_L_p,RTD_L_p)
Propdr2=c(LBuProp_adr2,RBuProp_adr2)
Propp=c(LBuProp_p,RBuProp_p)
Thetasdr2=c(LThetasFromPG_adr2,RThetasFromPG_adr2)
ThetasP=c(LThetasFromPG_p,RThetasFromPG_p)
# fdr each
BUp_f=p.adjust(BUp,method='fdr')
TDp_f=p.adjust(TDp,method='fdr')
Propp_f=p.adjust(Propp,method='fdr')
Thetasp_f=p.adjust(ThetasP,method='fdr')
# mask dr2s accordingly
BUdr2[BUp_f>0.05]=0
TDdr2[TDp_f>0.05]=0
Propdr2[Propp_f>0.05]=0
Thetasdr2[Thetasp_f>0.05]=0

# uncombine: seperate vecs for sep. hemis
BU_L=BUdr2[1:length(LBU_L_adr2)]
BU_R=BUdr2[(length(LBU_L_adr2)+1):length(BUdr2)]
TD_L=TDdr2[1:length(LBU_L_adr2)]
TD_R=TDdr2[(length(LBU_L_adr2)+1):length(BUdr2)]
Prop_L=Propdr2[1:length(LBU_L_adr2)]
Prop_R=Propdr2[(length(LBU_L_adr2)+1):length(BUdr2)]
Th_L=Thetasdr2[1:length(LBU_L_adr2)]
Th_R=Thetasdr2[(length(LBU_L_adr2)+1):length(BUdr2)]

# print out each for matlab friendly reading
write.table(BU_L,'~/results/PWs/FDRed_BU_L.csv',col.names=F,row.names=F,quote=F)
write.table(BU_R,'~/results/PWs/FDRed_BU_R.csv',col.names=F,row.names=F,quote=F)

write.table(TD_L,'~/results/PWs/FDRed_TD_L.csv',col.names=F,row.names=F,quote=F)
write.table(TD_R,'~/results/PWs/FDRed_TD_R.csv',col.names=F,row.names=F,quote=F)

write.table(Prop_L,'~/results/PWs/FDRed_Prop_L.csv',col.names=F,row.names=F,quote=F)
write.table(Prop_R,'~/results/PWs/FDRed_Prop_R.csv',col.names=F,row.names=F,quote=F)

write.table(Th_L,'~/results/PWs/FDRed_Th_L.csv',col.names=F,row.names=F,quote=F)
write.table(Th_R,'~/results/PWs/FDRed_Th_R.csv',col.names=F,row.names=F,quote=F)


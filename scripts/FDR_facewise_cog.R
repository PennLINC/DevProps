# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold dr2 values accordingly

# load in facewise output from left hemi
LBuProp_edr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_edr2.rds')
LBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_ep.rds')

# and right hemi
RBuProp_edr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_edr2.rds')
RBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_ep.rds')

# combine each
Propdr2=c(LBuProp_edr2,RBuProp_edr2)
Propp=c(LBuProp_p,RBuProp_p)
# fdr each
Propp_f=p.adjust(Propp,method='fdr')

# mask dr2s accordingly
Propdr2[Propp_f>0.05]=0

# uncombine: seperate vecs for sep. hemis
Prop_L=Propdr2[1:length(LBuProp_edr2)]
Prop_R=Propdr2[(length(LBuProp_edr2)+1):length(Propp_f)]

# print out each for matlab friendly reading
write.table(Prop_L,'~/results/PWs/FDRed_Prop_L_e.csv',col.names=T,row.names=F,quote=F)
write.table(Prop_R,'~/results/PWs/FDRed_Prop_R_e.csv',col.names=T,row.names=F,quote=F)


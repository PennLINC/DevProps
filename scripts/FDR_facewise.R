# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold dr2 values accordingly

# load in facewise output from left hemi
LBuProp_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_adr2.rds')
LBP_s=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_sdr2.rds')

LBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_p.rds')
LBP_sp=readRDS('/cbica/projects/pinesParcels/results/PWs/LBUProp_sp.rds')

# and right hemi
RBuProp_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_adr2.rds')
RBP_s=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_sdr2.rds')

RBuProp_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_p.rds')
RBP_sp=readRDS('/cbica/projects/pinesParcels/results/PWs/RBUProp_sp.rds')

# combine each
Propdr2=c(LBuProp_adr2,RBuProp_adr2)
Propp=c(LBuProp_p,RBuProp_p)
BUs=c(LBP_s,RBP_s[1:4842])
BUsP=c(LBP_sp,RBP_sp[1:4842])

# fdr each
Propp_f=p.adjust(Propp,method='fdr')
BUsP_f=p.adjust(BUsP,method='fdr')

# mask dr2s accordingly
Propdr2[Propp_f>0.05]=0
Propdr2[is.na(Propp_f)]=0

# uncombine: seperate vecs for sep. hemis
Prop_L=Propdr2[1:length(LBuProp_adr2)]
Prop_R=Propdr2[(length(LBuProp_adr2)+1):length(Propp)]

# print out each for matlab friendly reading
write.table(Prop_L,'~/results/PWs/FDRed_Prop_L.csv',col.names=F,row.names=F,quote=F)
write.table(Prop_R,'~/results/PWs/FDRed_Prop_R.csv',col.names=F,row.names=F,quote=F)


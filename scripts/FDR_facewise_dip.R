# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold dr2 values accordingly

# load in facewise output from left hemi
Ldip_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/LVerts_DR2.rds')
Ldip_p=readRDS('/cbica/projects/pinesParcels/results/PWs/LVerts_p.rds')

# and right hemi
Rdip_adr2=readRDS('/cbica/projects/pinesParcels/results/PWs/RVerts_DR2.rds')
Rdip_p=readRDS('/cbica/projects/pinesParcels/results/PWs/RVerts_p.rds')

# combine each
dr2s=c(Ldip_adr2,Rdip_adr2)
ps=c(Ldip_p,Rdip_p)

# fdr each
ps_f=p.adjust(ps,method='fdr')

# mask dr2s accordingly
dr2s[ps_f>0.05]=0

# uncombine: seperate vecs for sep. hemis
Dip_L=dr2s[1:length(Ldip_adr2)]
Dip_R=dr2s[(length(Ldip_adr2)+1):length(dr2s)]

# print out each for matlab friendly reading
write.table(Dip_L,'~/results/PWs/FDRed_Dip_L.csv',col.names=T,row.names=F,quote=F)
write.table(Dip_R,'~/results/PWs/FDRed_Dip_R.csv',col.names=T,row.names=F,quote=F)


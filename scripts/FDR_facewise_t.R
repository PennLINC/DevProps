# intake data from across hemispheres, fdr the p-distribution from both hemis, threshold t values accordingly

# load in facewise output from left hemi
L_t=readRDS('/cbica/projects/pinesParcels/results/PWs/r_vs_c_t_L.rds')
L_p=readRDS('/cbica/projects/pinesParcels/results/PWs/r_vs_c_p_L.rds')

# and right hemi
R_t=readRDS('/cbica/projects/pinesParcels/results/PWs/r_vs_c_t_R.rds')
R_p=readRDS('/cbica/projects/pinesParcels/results/PWs/r_vs_c_p_R.rds')

# combine each
ts=c(L_t,R_t)
ps=c(L_p,R_p)

# fdr each
ps_f=p.adjust(ps,method='fdr')

# mask dr2s accordingly
ts[ps_f>0.05]=0

# uncombine: seperate vecs for sep. hemis
ts_L=ts[1:length(L_t)]
ts_R=ts[(length(L_t)+1):length(ts)]

# print out each for matlab friendly reading
write.table(ts_L,'~/results/PWs/FDRed_ts_L.csv',col.names=F,row.names=F,quote=F)
write.table(ts_R,'~/results/PWs/FDRed_ts_R.csv',col.names=F,row.names=F,quote=F)

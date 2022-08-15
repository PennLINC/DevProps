#################################################################
# Calculate HierarchD, scatterplot Distances (FC vs. CFC, Yo vs. Old, Euc D vs. HierarchD)
#################################################################
#################################################################
###### load in Hierarchy, Calculate Distance, save to RDS for local loadin
vertPG_L<-read.csv('/cbica/projects/pinesParcels/results/PWs/vertPG_left.csv')
vertPG_R<-read.csv('/cbica/projects/pinesParcels/results/PWs/vertPG_right.csv')
# faces
facePG_L<-read.csv('/cbica/projects/pinesParcels/results/PWs/facePG_left.csv')
facePG_R<-read.csv('/cbica/projects/pinesParcels/results/PWs/facePG_right.csv')
# Initialize HierarchD matrices
HD_v_l<-matrix(nrow=dim(vertPG_L)[1],ncol=dim(vertPG_L)[1])
HD_v_r<-matrix(nrow=dim(vertPG_R)[1],ncol=dim(vertPG_R)[1])
HD_f_l<-matrix(nrow=dim(facePG_L)[1],ncol=dim(facePG_L)[1])
HD_f_r<-matrix(nrow=dim(facePG_R)[1],ncol=dim(facePG_R)[1])
# medial wall vertices
mw_l_verts=as.logical(read.csv('~/data/mw_boolean_l.csv'))
mw_r_verts=as.logical(read.csv('~/data/mw_boolean_r.csv'))
##### Calculate HierarchD

# vertices

# L
for (i in 1:dim(vertPG_L)[1]){
  for (j in 1:dim(vertPG_L)[1]){
    # absolute valued distance
    HD_v_l[i,j]<-abs(vertPG_L[i,1]-vertPG_L[j,1])
  }
}

# mask medial wall
HD_v_l=HD_v_l[mw_l_verts,mw_l_verts]
# upper tri is only hope of shrinking file enough to load into mount for plotting
HD_v_l=HD_v_l[upper.tri(HD_v_l)]
# saveout
saveRDS(HD_v_l,'/cbica/projects/pinesParcels/data/HD_v_l.rds')

# R

for (i in 1:dim(vertPG_R)[1]){
  print('right')
  for (j in 1:dim(vertPG_R)[1]){
    # absolute valued distance
    HD_v_r[i,j]<-abs(vertPG_R[i,1]-vertPG_R[j,1])
  }
}

print('done with vertices')

# mask medial wall
HD_v_r=HD_v_r[mw_r_verts,mw_r_verts]
# upper tri is only hope of shrinking file enough to load into mount for plotting
HD_v_r=HD_v_r[upper.tri(HD_v_r)]
# saveout
HD_v_r=HD_v_r[upper.tri(HD_v_r)]
saveRDS(HD_v_r,'/cbica/projects/pinesParcels/data/HD_v_r.rds')

# faces

# L
for (i in 1:dim(facePG_L)[1]){
  for (j in 1:dim(facePG_L)[1]){
    # absolute valued distance
    HD_f_l[i,j]<-abs(facePG_L[i,1]-facePG_L[j,1])
  }
}

HD_f_l=HD_f_l[upper.tri(HD_f_l)]
# saveout
saveRDS(HD_f_l,'/cbica/projects/pinesParcels/data/HD_f_l.rds')

# R 

for (i in 1:dim(facePG_R)[1]){
  print(i)
  for (j in 1:dim(facePG_R)[1]){
    # absolute valued distance
    HD_f_r[i,j]<-abs(facePG_R[i,1]-facePG_R[j,1])
  }
}

HD_f_r=HD_f_r[upper.tri(HD_f_r)]
# saveout
saveRDS(HD_f_r,'/cbica/projects/pinesParcels/data/HD_f_r.rds')

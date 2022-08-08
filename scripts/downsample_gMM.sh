# 32k PG (jk its a myelin map)
PG32=~/data/HCD628_Winter2021.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii

# separate hemispheres - left
wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_LEFT ~/data/hcpd.mm_L.func.gii

# right hemi
wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_RIGHT ~/data/hcpd.mm_R.func.gii

### resample both hemis to 10k vertices
# left hemisphere
wb_command -metric-resample ~/data/hcpd.mm_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ~/data/hcpd.mm_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample ~/data/hcpd.mm_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ~/data/hcpd.mm_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii



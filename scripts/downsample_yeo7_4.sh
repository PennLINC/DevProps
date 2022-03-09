# 32k PG
PG32=/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii

# separate hemispheres - left
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_LEFT /cbica/projects/abcdfnets/data/hcp.gradients_L.func.gii 

# right hemi
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_RIGHT /cbica/projects/abcdfnets/data/hcp.gradients_R.func.gii

### resample both hemis to 3k vertices
# left hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/abcdfnets/data/hcp.gradients_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA /cbica/projects/abcdfnets/data/hcp.gradients_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/abcdfnets/data/hcp.gradients_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA /cbica/projects/abcdfnets/data/hcp.gradients_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii


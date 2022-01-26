# 32k PG
PG32=/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii

# separate hemispheres - left
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_LEFT /cbica/projects/abcdfnets/data/hcp.gradients_L.func.gii 

# right hemi
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_RIGHT /cbica/projects/abcdfnets/data/hcp.gradients_R.func.gii

### resample both hemis to 1k vertices - BARY centric! 
# left hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/abcdfnets/data/hcp.gradients_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ~/dropbox/lh.sphere3.surf.gii BARYCENTRIC /cbica/projects/abcdfnets/data/hcp.gradients_L_1k.func.gii

# right hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/abcdfnets/data/hcp.gradients_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ~/dropbox/rh.sphere3.surf.gii BARYCENTRIC /cbica/projects/abcdfnets/data/hcp.gradients_R_1k.func.gii

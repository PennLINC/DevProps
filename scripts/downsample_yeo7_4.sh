# yeo7
Y7=/cbica/projects/pinesParcels/data/Yeo7_Yeo-Buckner-Choi-Raut.dlabel.nii

# separate hemispheres - left
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $Y7 COLUMN -label CORTEX_LEFT /cbica/projects/pinesParcels/data/y7_L.label.gii 

# right hemi
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $Y7 COLUMN -label CORTEX_RIGHT /cbica/projects/pinesParcels/data/y7_R.label.gii

### resample both hemis to 3k vertices
# left hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -label-resample /cbica/projects/pinesParcels/data/y7_L.label.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/data/y7_L_3k.label.gii -area-metrics ~/dropbox/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -label-resample /cbica/projects/pinesParcels/data/y7_R.label.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/data/y7_R_3k.label.gii -area-metrics ~/dropbox/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii ~/dropbox/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii


# subject name
subj=$1

# output filepath
OFP=/cbica/projects/pinesParcels/results/PWs/Proced/${subj}/

# subject's principal gradients
LPGfp=${OFP}${subj}_PG_L_10k_rest.func.gii
RPGfp=${OFP}${subj}_PG_R_10k_rest.func.gii

# left hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample ${LPGfp} /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA ${OFP}${subj}_PG_L_32k.func.gii -area-metrics /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii
# right hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample ${RPGfp} /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA ${OFP}${subj}_PG_R_32k.func.gii -area-metrics /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii

# combine
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-create-dense-from-template /cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients.dscalar.nii ${OFP}${subj}_PG_LR_32k_rest.dscalar.nii -metric CORTEX_LEFT ${OFP}${subj}_PG_L_32k.func.gii -metric CORTEX_RIGHT ${OFP}${subj}_PG_R_32k.func.gii

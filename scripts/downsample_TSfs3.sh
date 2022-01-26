# subject name
subj=$1

# subject's aggregated time series
AgTS=/cbica/projects/hcpd/data/motMasked_contSegs/${subj}/${subj}_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii

# separate hemispheres - left
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_L_AggTS.func.gii 

# right hemi
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_R_AggTS.func.gii

### resample both hemis to 10k vertices
# left hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_L_AggTS.func.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/lh.sphere3.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_L_1k.func.gii -area-metrics /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/hcp.gradients_L_1k.func.gii

# right hemisphere
/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_R_AggTS.func.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/rh.sphere3.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_R_1k.func.gii -area-metrics /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/standard_mesh_atlases/resample_fsaverage/hcp.gradients_R_1k.func.gii

# convert to mgh for reading individual hemisphere time series into matlab
mri_convert /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_L_1k.func.gii /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_L_1k.mgh

mri_convert /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_R_1k.func.gii /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_R_1k.mgh

# cp to hcpd

cp /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_L_1k.mgh /cbica/projects/hcpd/data/motMasked_contSegs/${subj}/

cp /cbica/projects/pinesParcels/results/PWs/PreProc/${subj}/${subj}_AggTS_R_1k.mgh /cbica/projects/hcpd/data/motMasked_contSegs/${subj}/

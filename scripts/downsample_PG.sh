subj=$1

# 32k PG
#PG32=/cbica/projects/pinesParcels/results/PWs/Proced/${subj}/${subj}_PG_LR_32k_rest.dscalar.nii 

# separate hemispheres - left
#wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_LEFT ~/tmp/gradients_L.func.gii 

# right hemi
#wb_command -cifti-separate $PG32 COLUMN -metric CORTEX_RIGHT ~/tmp/gradients_R.func.gii

### resample both hemis to 3k vertices
# left hemisphere
wb_command -metric-resample /cbica/projects/pinesParcels/results/PWs/Proced/${subj}/${subj}_PG_L_10k_rest.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/results/PWs/Proced/${subj}/${subj}_gradients_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample /cbica/projects/pinesParcels/results/PWs/Proced/${subj}/${subj}_PG_R_10k_rest.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA /cbica/projects/pinesParcels/results/PWs/Proced/${subj}/${subj}_gradients_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii


# subject name
subj=$1

# subject's aggregated time series
AgTS=/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_AggTS.dtseries.nii

# separate hemispheres - left
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_L_AggTS.func.gii 

# right hemi
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_R_AggTS.func.gii

### resample both hemis to 10k vertices
# left hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_L_AggTS.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_AggTS_L_10k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii

# right hemisphere
/cbica/projects/abcdfnets/scripts/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_R_AggTS.func.gii ~/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_AggTS_R_10k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

# written out as resamp_32K_network${i}_L.32K_fs_LR.func.gii


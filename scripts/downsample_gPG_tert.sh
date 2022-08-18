### resample both hemis to 3k vertices

### young
# left hemisphere
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_young_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ~/data/hcpd.youngPG_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_young_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ~/data/hcpd.youngPG_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

### mid
# left hemisphere
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_mid_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ~/data/hcpd.midPG_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_mid_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ~/data/hcpd.midPG_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

### old
# left hemisphere
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_old_L.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ~/data/hcpd.oldPG_L_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right
wb_command -metric-resample /cbica/projects/pinesParcels/data/princ_gradients/hcpd_old_R.func.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ~/data/hcpd.oldPG_R_3k.func.gii -area-metrics ~/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii ~/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii


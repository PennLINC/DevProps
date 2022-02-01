# NOTE: THIS MUST BE RUN BY HCPD USER
# load in subj list, loop over
cat /cbica/projects/hcpd/dropbox/G600TRs.txt | while read line; do 

# just to remove the newline
Subj=$(echo $line | cut -b -11);

# for each network i at k=17
for i in {1..17}; do

# seperate the hemis
/cbica/projects/hcpd/dropbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/AtlasLoading_Network_${i}.dscalar.nii COLUMN -metric CORTEX_LEFT /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/lh_Network_${i}.func.gii 

/cbica/projects/hcpd/dropbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -cifti-separate /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/AtlasLoading_Network_${i}.dscalar.nii COLUMN -metric CORTEX_RIGHT /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/rh_Network_${i}.func.gii

# resample soft parcel to hcp-surface
# left hemisphere
/cbica/projects/hcpd/dropbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/lh_Network_${i}.func.gii /cbica/projects/hcpd/dropbox/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /cbica/projects/hcpd/dropbox/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/resamp_3k_network_L${i}.func.gii -area-metrics /cbica/projects/hcpd/dropbox/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /cbica/projects/hcpd/dropbox/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
/cbica/projects/hcpd/dropbox/workbench/workbench-1.2.3/exe_rh_linux64/wb_command -metric-resample /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/rh_Network_${i}.func.gii /cbica/projects/hcpd/dropbox/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /cbica/projects/hcpd/dropbox/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA /cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/${Subj}/resamp_3k_network_R${i}.func.gii -area-metrics /cbica/projects/hcpd/dropbox/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /cbica/projects/hcpd/dropbox/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# end over networks loop
done
# end over subject loop
done

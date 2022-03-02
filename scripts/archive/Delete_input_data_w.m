function Delete_input_data(subj)
sname=subj;
DL_dir=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/'];
deleteCommand=['rm -r ' DL_dir ' -f'];
system(deleteCommand) 
subjTxtFile=['/cbica/projects/abcdfnets/nda-abcd-s3-downloader/' subj '.txt'];
deleteCommand=['rm ' subjTxtFile ' -f'];
system(deleteCommand) 


% now delete NMF output for extra lean folders
sbjDataFile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/sbjData.mat'];
deleteCommand=['rm ' sbjDataFile ' -f'];
system(deleteCommand)
sbjIterRes=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/Iteration_res.mat'];
deleteCommand=['rm ' sbjIterRes ' -f'];
system(deleteCommand)
sbjPreProc=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/Preprocessing.mat'];
deleteCommand=['rm ' sbjPreProc ' -f'];
system(deleteCommand)
initUVFile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/init_UV.mat'];
deleteCommand=['rm ' initUVFile ' -f'];
system(deleteCommand)
CentroidFile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/res_cen.mat'];
deleteCommand=['rm ' CentroidFile ' -f'];
system(deleteCommand)



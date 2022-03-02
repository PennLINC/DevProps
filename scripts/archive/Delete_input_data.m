function Delete_input_data(subj)
sname=subj;
DL_dir=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/'];
deleteCommand=['rm -r ' DL_dir ' -f'];
system(deleteCommand) 
subjTxtFile=['/cbica/projects/abcdfnets/nda-abcd-s3-downloader/' subj '.txt'];
deleteCommand=['rm ' subjTxtFile ' -f'];
system(deleteCommand) 

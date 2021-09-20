import numpy as np
import nibabel as nb
import sys
import nilearn
subj = sys.argv[1]
# temp file dir
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# initialize 25 x 59142 np
dataArray=np.empty([25,91282])
for f in range(25):
	# read in combined LR frame
	subjDatafn=parentfp+str(subj)+'_PG_LR_f'+str(f)+'_32k.dscalar.nii'
	subjData=nb.load(subjDatafn)
	# insert into dataframe
	dataArray[f,:]=subjData.dataobj[0,:]

# create cifti from dataframe and templates
timeAxis=nb.cifti2.SeriesAxis(start=0,size=25,step=1)
spaceAxis=subjData.header.get_axis(1)
new_img=nb.Cifti2Image(dataArray,(timeAxis,spaceAxis))
filename='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/' + str(subj) + '_QPP.dtseries.nii'
new_img.to_filename(filename)

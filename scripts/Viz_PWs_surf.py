import nibabel as nb
import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
import subprocess
# specify subject
subj = sys.argv[1]
# set parent fp
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
wavedir='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# set plot output fp
pofp='/cbica/projects/abcdfnets/results/PWplots/' + str(subj)
os.mkdir(pofp)
###### import TS
TSfp=parentfp + str(subj) + '_downsamp_rest.npy'
TS=np.load(TSfp)
# normalize
TS=TS-np.mean(TS,0)
STDS=np.std(TS)
TS=TS/STDS
###### import PW instances
PWinstfp=wavedir + str(subj) + '_PW1_peaks.npy'
PWinstances=np.load(PWinstfp)
# for each instance, plot its frames
for i in range(len(PWinstances)):
	# get startframe
	startFrame=PWinstances[i]
	EndFrame=startFrame+25
	# extract from TS
	PW=TS[:,int(startFrame):int(EndFrame)]
	### for each frame, print out
	for frame in range(25):
		# save wave time series frame as csv for matlab to read
		wffn=parentfp+str(subj)+'waveTS.csv'
		wf=TS[:,int(startFrame+frame)]
		np.savetxt(wffn,wf,delimiter=',')
		outputName=pofp+'/'+str(subj)+'_PWinst_'+str(i)+'_frame_'+str(frame)
		# call matlab script
		matlabCommand = 'PBP_vertWiseEffect(' + "'" + wffn + "','" + outputName + "')"		
		#bashCommand = "matlab -nodisplay -r " + matlabCommand
		#print(bashCommand)
		subprocess.run(["matlab", "-nodisplay","-r",matlabCommand])	

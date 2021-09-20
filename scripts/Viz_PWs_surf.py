import nibabel as nb
import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
# specify subject
subj = sys.argv[1]
# set parent fp
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
wavedir='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# set plot output fp
pofp='/cbica/projects/abcdfnets/results/PWplots/' + str(subj)
os.mkdir(pofp)
###### import TS
TSfp=parentfp + str(subj) + '_p2mm_masked_filtered_rest.dtseries.nii'
TS=nb.load(TSfp)
# extract hemis
CL=TS.dataobj[:,hcp.struct.cortex_left]
# each vertex normalized w/r/t global SD
# L
Lstds=np.std(CL)
CL=CL/Lstds;
###### import PW instances
PWinstfp=wavedir + str(subj) + '_PW1_peaks.npy'
PWinstances=np.load(PWinstfp)
# for each instance, plot its frames
for PWinst in range(len(PWinst)):
	# get startframe
	startFrame=PWinstances[PWinst]
	EndFrame=startFrame+25
	# extract from TS
	PW_L=CL[int(startFrame):int(EndFrame),:]
	### for each frame, print out
	for frame in range(25):
		display=plotting.view_surf(hcp.mesh.inflated_left, hcp.left_cortex_data(CL[frame,:]), threshold=0, bg_map=hcp.mesh.sulc_left)

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
###### import GS
GS=np.genfromtxt(parentfp + str(subj) '_p2mm_masked_filtered_rest_GS.csv',delimiter=',')
###### import PG
PGfp=wavedir + str(subj) + '_PG_LR_32k.dscalar.nii'
PG=nb.load(PGfp)
###### import TS
TSfp=parentfp + str(subj) + '_p2mm_masked_filtered_rest.dtseries.nii'
TS=nb.load(TSfp)
# extract hemis
CL=TS.dataobj[:,hcp.struct.cortex_left]
CR=TS.dataobj[:,hcp.struct.cortex_right]
# each vertex normalized w/r/t global SD
# L
Lstds=np.std(CL)
CL=CL/Lstds;
CLt=np.transpose(CL)
# R
Rstds=np.std(CR)
CR=CR/Rstds;
CRt=np.transpose(CR)
###### import PW instances
PWinstfp=wavedir + str(subj) + '_PW1_peaks.npy'
PWinst=np.load(PWinstfp)
###### import PW 
PWfp
###### plot
###### save plot to subj folder

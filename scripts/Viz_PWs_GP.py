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
PGL=pg.dataobj[0,hcp.struct.cortex_left]
# and right 1st pg
PGR=pg.dataobj[0,hcp.struct.cortex_right]
# to organize vertices in terms of their position on the PG
PGindicesL=np.argsort(PGL)
PGindicesR=np.argsort(PGR)
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
# R
Rstds=np.std(CR)
CR=CR/Rstds;
# sort cortical data /w/r/t pg
CL=CL[:,PGindicesL]
CR=CR[:,PGindicesR]
# transpose for grayplots
CLt=np.transpose(CL)
CRt=np.transpose(CR)
###### import PW instances
PWinstfp=wavedir + str(subj) + '_PW1_peaks.npy'
PWinst=np.load(PWinstfp)
###### plot
# calculate data volume:100 aspect ratio correction number
ARCN=(CL.shape[0])/100
# Left Hemi
plt.subplot(3,1,2)
for x in range(len(PWinst)):
	plt.axvline(x=PWinst[x]/ARCN,linewidth=.2,dashes=[1,3,1,3])
plt.matshow(CLt,cmap='gray',vmin=-2, vmax=2,extent=[0,100,0,1], aspect=15.5,fignum=False)

# Right Hemi
plt.subplot(3,1,3)
for x in range(len(PWinst)):
        plt.axvline(x=PWinst[x]/ARCN,linewidth=.2,dashes=[1,3,1,3])
plt.matshow(CRt,cmap='gray',vmin=-2, vmax=2,extent=[0,100,0,1], aspect=15.5,fignum=False)

# plot GS
plt.subplot(3,1,1)
for x in range(len(PWinst)):
        plt.axvline(x=PWinst[x],linewidth=.2,dashes=[1,3,1,3])
plt.plot(GS,c='black')
###### save plot to subj folder
plt.savfig(pofp + '_PWs_GP.png',bbox_inches='tight')

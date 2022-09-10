import nibabel as nb
import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
# output fp
pofp = '/cbica/projects/pinesParcels/PWs/scripts/'
# prevent plt from trying to render mid script, only saveout
matplotlib.use('agg')
# set resolution
plt.rcParams['figure.dpi'] = 700
# set parent fp
parentfp='/cbica/projects/pinesParcels/results/PWs/PreProc/' + str(subj) + '/'
### load in lobe indices
leftLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_lh.csv')
rightLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_rh.csv')
###### import PG
PGfp='/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'
pgL=nb.load(PGfp)
PGfp='/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'
pgR=nb.load(PGfp)
# to organize vertices in terms of their position on the PG
PGindicesL=np.argsort(pgL.agg_data()[0])
PGindicesR=np.argsort(pgR.agg_data()[0])
###### import TS
TSfp=parentfp + str(subj) + '_AggTS_L_3k.func.gii'
gif_CL=nb.load(TSfp)
CL=np.array(gif_CL.agg_data()[:])
TSfp=parentfp + str(subj) + '_AggTS_R_3k.func.gii'
gif_CR=nb.load(TSfp)
CR=np.array(gif_CR.agg_data()[:])
#### extract time series and PG indices for frontal lobe
# get indices in appropriate format
lobeCode=6500
lobeIndexL=np.where(leftLubs==lobeCode)
lobeIndexR=np.where(rightLubs==lobeCode)
# extract
CL=CL[:,lobeIndexL]
CR=CR[:,lobeIndexR]
PGindicesL=PGindicesL[lobeIndexL]
PGindicesR=PGindicesR[lobeIndexR]
CL=np.squeeze(CL)
CR=np.squeeze(CR)
# re-do ordering for within lobe
PGindicesL=np.argsort(PGindicesL)
PGindicesR=np.argsort(PGindicesR)
# each vertex normalized w/r/t global SD
# L
Lstds=np.std(CL)
CL=CL/Lstds;
# R
Rstds=np.std(CR)
CR=CR/Rstds;
# sort cortical data /w/r/t pg
CL=CL[800:1300,PGindicesL]
CR=CR[800:1300,PGindicesR]
# transpose for grayplots
CLt=np.transpose(CL)
CRt=np.transpose(CR)

###### plot
# calculate data volume:100 aspect ratio correction number
ARCN=(CL.shape[0])/100
# Left Hemi
plt.subplot(2,1,1)
plt.matshow(CLt,cmap='gray',vmin=-1.5, vmax=1.5,extent=[0,100,0,1], aspect=50,fignum=False)
# Right Hemi
plt.subplot(2,1,2)
plt.matshow(CRt,cmap='gray',vmin=-1.5, vmax=1.5,extent=[0,100,0,1], aspect=50,fignum=False)
###### save plot to subj folder
plt.savefig(pofp + 'test_GPtall.png',bbox_inches='tight')

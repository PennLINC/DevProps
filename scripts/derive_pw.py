import scipy
import sys
import hcp_utils as hcp
import numpy as np
import nibabel as nb
import csv
import tables
from scipy.io import loadmat
import nibabel as nb
import numpy as np
from QPP_cpac_edit import detect_qpp
subj = sys.argv[1]
# parent filepath for time series
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# load in downsampled aggregate time series
AgTS_l=nb.load(parentfp+subj+'_AggTS_L_10k.func.gii')
# data array
LdArr=AgTS_l.darrays
# right
AgTS_r=nb.load(parentfp+subj+'_AggTS_R_10k.func.gii')
# data array
RdArr=AgTS_r.darrays
# initialize array to insert time series into, fsaverage5 resolution
timeSeries=np.zeros((20484,len(LdArr)))
# populate array with both hemispheres
for TR in range(len(LdArr)):
	# odd error where it makes me broadcast into first 10243 spots (0 inclusive), overwrite with second command
	timeSeries[0:10242,TR]=LdArr[TR].data
	timeSeries[10242:20484,TR]=RdArr[TR].data
# find where average and SD over all TRs is 0
Averages=np.average(timeSeries,axis=1)
StDs=np.std(timeSeries,axis=1)
### remove medial wall
MWindices = np.where((Averages == 0) | (StDs == 0))
MWindices=MWindices[0]
NoMWts=np.delete(timeSeries,MWindices,axis=0)
nonMWindices = np.where((Averages != 0) | (StDs != 0))
nonMWindices = nonMWindices[0]
# get valid resting segment for now
segmentfnTr=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'
segmentTr=np.genfromtxt(segmentfnTr,delimiter=',')
# calc # of resting TRs
lastStart=segmentTr[-1,0]
lestLength=segmentTr[-1,1]
RSlength=lastStart+lestLength-1
# extract resting segment
restSeg=NoMWts[:,0:int(RSlength)]
# save out this segment for PW viz later
# run QPP
QPPout=detect_qpp(segmentTr,restSeg,1,25,10,.1,10,convergence_iterations=10)
# save 1st qpp, metrics, and instances
subjPWfn_peaks=childfp+str(subj)+'_PW1_peaks'
np.save(subjPWfn_peaks,QPPout[1])
subjPWfn_metrics=childfp+str(subj)+'_PW1_metrics'
np.save(subjPWfn_metrics,QPPout[2])
### convert to cifti format for upsampling
# MW included array for saveout
subjPWMW=np.zeros((20484,25))
subjPWMW[nonMWindices,:]=QPPout[0]
# and MW included for greater time series saveout as well
restSegOut=np.zeros((20484,restSeg.shape[1]))
restSegOut[nonMWindices,:]=restSeg
np.save(parentfp+str(subj)+'_downsamp_rest',restSegOut)
# new cifti axis, 25 to match QPP
#newAxis=nb.cifti2.SeriesAxis(start=0,size=25,step=1)
# extract spatial axis from downsampled template
#GradsL=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_L_10k.func.gii')
#GradsR=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_R_10k.func.gii')
# saveout each frame individually
#for f in range(25):
#	PWL=GradsL.darrays[0]
#	PWR=GradsR.darrays[0]
#	PWL.data=subjPWMW[0:10242,f]
#	PWR.data=subjPWMW[10242:20484,f]
#	GradsL.darrays[0]=PWL
#	GradsR.darrays[0]=PWR
#	Lfp=parentfp+str(subj)+'_QPP1_f'+str(f)+'_L_10k.func.gii'
#	Rfp=parentfp+str(subj)+'_QPP1_f'+str(f)+'_R_10k.func.gii'
#	nb.save(GradsL,Lfp)
#	nb.save(GradsR,Rfp)

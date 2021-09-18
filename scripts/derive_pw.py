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
from QPP_cpac import detect_qpp
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
### remove medial wall
# find where average and SD over all TRs is 0
Averages=np.average(timeSeries,axis=1)
StDs=np.std(timeSeries,axis=1)
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
# run QPP
QPPout=detect_qpp(NoMWts,1,25,1,.1,1,convergence_iterations=1)
# save 1st qpp, metrics, and instances
subjPWfn_verts=childfp+str(subj)+'_PW1_verts'
subjPWfn_peaks=childfp+str(subj)+'_PW1_peaks'
subjPWfn_metrics=childfp+str(subj)+'_PW1_metrics'
### convert to cifti format for upsampling
# MW included array for saveout
subjPWMW=np.zeros((20484,25))
subjPWMW[nonMWindices,:]=QPPout[0]
# new cifti axis, 25 to match QPP
newAxis=nb.cifti2.SeriesAxis(start=0,size=25,step=1)
# loading in downsampled PG as template
GradsL=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_L_10k.func.gii')
GradsR=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_R_10k.func.gii')
# extract spatial axis from downsampled template
LgOrdAx=GradsL.header.get_axis(1)
RgOrdAx=GradsR.header.get_axis(1)
# extract left and right hemi from medial wall-filled object, implant header info from templates
new_imgL=nb.Cifti2Image(subjPWMW[0:10242,:],(newAxis,LgOrdAx))
new_imgR=nb.Cifti2Image(subjPWMW[10242:20484,:],(newAxis,RgOrdAx))
nb.save(new_imgL,subjPWfn_verts+'_L.func.gii')
nb.save(new_imgR,subjPWfn_verts+'_R.func.gii')
np.save(subjPWfn_verts,QPPout[0])
np.save(subjPWfn_peaks,QPPout[1])
np.save(subjPWfn_metrics,QPPout[2])

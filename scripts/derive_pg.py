import scipy
import sys
import hcp_utils as hcp
import numpy as np
import nibabel as nb
import csv
import tables
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances
from mapalign import embed
from sklearn.metrics.pairwise import cosine_similarity

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
# create FC matrix
fcmatrix=np.corrcoef(NoMWts)
# Get number of nodes
N = fcmatrix.shape[0]
# Generate percentile thresholds for 90th percentile
perc = np.array([np.percentile(x, 90) for x in fcmatrix])
# Threshold each row of the matrix by setting values below 90th percentile to 0
for i in range(fcmatrix.shape[0]):
	fcmatrix[i, fcmatrix[i,:] < perc[i]] = 0

# Check for minimum value
print("Minimum value is %f" % fcmatrix.min())
# set negative values to 0
fcmatrix[fcmatrix < 0] = 0
# Now we are dealing with sparse vectors. Cosine similarity is used as affinity metric
#aff = cosine_similarity_n_space(fcmatrix)
aff = cosine_similarity(fcmatrix)
# now compute diffusion map
emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)
# select PG based on max absolute corr with group-level PG
# (loading in downsampled PG)
GradsL=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_L_10k.func.gii')
GradsR=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_R_10k.func.gii')
# extract 1st grad
PGL=GradsL.darrays[0]
PGR=GradsR.darrays[0]
# array for pg1
PG1=np.zeros((20484,1))
PG1[0:10242,0]=PGL.data
PG1[10242:20484,0]=PGR.data
# same MW mask 
NoMWPG=np.delete(PG1,MWindices,axis=0)
# initialize correlation list for each of this subject's first five gradients
gradCors=np.zeros(5)
# checking first five gradients is probably overkill but there are a lot of abcd subjs to be potential edge cases
for g in range(5):
	spear_man=scipy.stats.spearmanr(NoMWPG,emb[:,g])
	gradCors[g]=np.absolute(spear_man.correlation)

print(max(gradCors))
# argmax to find location of best PG-equivalent estimate
subjPG_ind=np.argmax(gradCors)
subjPG=emb[:,subjPG_ind]
# match gradient directionality (i.e., is transmodal pos and unimodal neg?)
if scipy.stats.spearmanr(subjPG,NoMWPG).correlation < 0:
	subjPG=subjPG * -1

# replace group PG values with this subj's for printout
# MW included array for saveout
subjPGMW=np.zeros((20484,1))
subjPGMW[nonMWindices,0]=subjPG
PGL.data=subjPGMW[0:10242,0]
PGR.data=subjPGMW[10242:20484,0]
GradsL.darrays[0]=PGL
GradsR.darrays[0]=PGR
# save out subject-specific PG for upsampling
LPGfp=parentfp+subj+'_PG_L_10k.func.gii'
RPGfp=parentfp+subj+'_PG_R_10k.func.gii'
nb.save(GradsL,LPGfp)
nb.save(GradsR,RPGfp)
# extract variance explained by this grad
lambdas=res['lambdas']/sum(res['lambdas'])
PGlambd=lambdas[subjPG_ind]
# save index of "PG" for later - which gradient matched THE principal gradient for this subject
# also save variance explained by this gradient
subjPGfn_index=childfp+str(subj)+'_PG1_index'
np.save(subjPGfn_index,[subjPG_ind,PGlambd])

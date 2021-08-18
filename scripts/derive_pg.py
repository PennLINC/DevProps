import sys
import numpy as np
import nibabel as nib
import csv
import tables
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances
from mapalign import embed

subj = sys.argv[1]

# parent filepath for time series
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'

# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'

# all the scan types
tasks=['rest','SST','nback','mid']

for T in range(len(tasks)):
	# grayordinate fc matrix (cortical only)
	Gfp=parentfp+str(subj)+'_FullFCmat_'+tasks[T]+'.npy'
	fcmatrix = np.load(Gfp)
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
	aff = 1 - pairwise_distances(fcmatrix, metric = 'cosine')
	# now compute diffusion map
	emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)
	PGfp=childfp+str(subj)+'_PGs_'+tasks[T]
	resFP=childfp+str(subj)+'_res_'+tasks[T]
	np.save(emb,PGfp)
	np.save(res,resFP)

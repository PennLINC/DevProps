import scipy
import sys
import numpy as np
import nibabel as nb
import csv
import tables
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances
from mapalign import embed
from sklearn.metrics.pairwise import cosine_similarity
from scipy.io import savemat

# load in aggregated matrices
outFN='/cbica/projects/pinesParcels/results/PWs/Gro_CFC_L.mat'
fcmatrixL=loadmat(outFN)
fcmatrixL=fcmatrixL['groL']
outFN='/cbica/projects/pinesParcels/results/PWs/Gro_CFC_R.mat'
fcmatrixR=loadmat(outFN)
fcmatrixR=fcmatrixR['groR']
#### LEFT
# Get number of nodes
N = fcmatrixL.shape[0]
# Generate percentile thresholds for 90th percentile
perc = np.array([np.percentile(x, 90) for x in fcmatrixL])
# Threshold each row of the matrix by setting values below 90th percentile to 0
for i in range(fcmatrixL.shape[0]):
	fcmatrixL[i, fcmatrixL[i,:] < perc[i]] = 0

# Check for minimum value
print("Minimum value is %f" % fcmatrixL.min())
# set negative values to 0
fcmatrixL[fcmatrixL < 0] = 0
# Now we are dealing with sparse vectors. Cosine similarity is used as affinity metric
#aff = cosine_similarity_n_space(fcmatrix)
aff = cosine_similarity(fcmatrixL)
# now compute diffusion map
emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)
savemat('/cbica/projects/pinesParcels/results/PWs/embL.mat',{'emb':emb})
### RIGHT
# Get number of nodes
N = fcmatrixR.shape[0]
# Generate percentile thresholds for 90th percentile
perc = np.array([np.percentile(x, 90) for x in fcmatrixR])
# Threshold each row of the matrix by setting values below 90th percentile to 0
for i in range(fcmatrixR.shape[0]):
        fcmatrixR[i, fcmatrixR[i,:] < perc[i]] = 0

# Check for minimum value
print("Minimum value is %f" % fcmatrixR.min())
# set negative values to 0
fcmatrixR[fcmatrixR < 0] = 0
# Now we are dealing with sparse vectors. Cosine similarity is used as affinity metric
#aff = cosine_similarity_n_space(fcmatrix)
aff = cosine_similarity(fcmatrixR)
# now compute diffusion map
emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)
savemat('/cbica/projects/pinesParcels/results/PWs/embR.mat',{'emb':emb})



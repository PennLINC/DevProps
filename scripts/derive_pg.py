import sys
import hcp_utils as hcp
import numpy as np
import nibabel as nib
import csv
import tables
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances
from mapalign import embed
from sklearn.metrics.pairwise import cosine_similarity

# define a function to chunk cosine affin. into smaller matrices, gets stuck otherwise
# https://stackoverflow.com/questions/40900608/cosine-similarity-on-large-sparse-matrix-with-numpy

def cosine_similarity_n_space(m, batch_size=200):
    ret = np.ndarray((m.shape[0], m.shape[0]))
    for row_i in range(0, int(m.shape[0] / batch_size) + 1):
        print(row_i)
	start = row_i * batch_size
        end = min([(row_i + 1) * batch_size, m.shape[0]])
        if end <= start:
            break
        rows = m[start: end]
        sim = cosine_similarity(rows, m) # rows is O(1) size
        ret[start: end] = sim
    return ret

subj = sys.argv[1]

# parent filepath for time series
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'

# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'

# grayordinate fc matrix (cortical only)
Gfp=parentfp+str(subj)+'_FullFCmat.npy'
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
aff = cosine_similarity_n_space(fcmatrix)
# now compute diffusion map
emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)
PGfp=childfp+str(subj)+'_PGs'
resFP=childfp+str(subj)+'_res'
np.save(emb,PGfp)
np.save(res,resFP)
# select PG based on max absolute corr with group-level PG
Grads=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii')
GradsCort=Grads.dataobj[:,hcp.struct.cortex]
gPG=GradsCort[0,:]
# initialize correlation list for each of this subject's first five gradients
gradCors=np.zeros(5)
# checking first five gradients is probably overkill but there are a lot of abcd subjs to be potential edge cases
for g in range(5):
	gradCors[g]=np.absolute(np.corrcoef(gPG,emb[:,g]))

print(max(gradCors))
# argmax to find location of best PG-equivalent estimate
subjPG_ind=np.argmax(gradCors)
subjPG=emb[:,subjPG_ind]
subjPGfn=childfp+str(subj)+'_PG1'
np.save(subjPG,subjPGfn)

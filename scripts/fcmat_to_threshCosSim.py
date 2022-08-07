#!/cbica/projects/pinesParcels/miniconda3/envs/difembs/bin/python
# -*- coding: utf-8 -*-
# lotta code taken from https://github.com/NeuroanatomyAndConnectivity/gradient_analysis
import sys
import numpy as np
import nibabel as nib
import csv
import tables
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances

# call with "s number" as only argument

snum = int(sys.argv[1])

# Load in subjects for filepaths
s_file = open("/cbica/projects/pinesParcels/data/bblids.txt")
sfile_contents = s_file. read()
scontents_split = sfile_contents. splitlines()
# get this specific subject
sid=scontents_split[snum]
# filepath of fcmat is
fcfp = "/cbica/projects/pinesParcels/data/CombinedData/" + str(sid) + "/vertexwise_fc_mat.csv"
# load in the mat file
fcmatrix= np.genfromtxt(fcfp,delimiter=',')

# Get number of nodes
N = fcmatrix.shape[0]

# Generate percentile thresholds for 90th percentile
perc = np.array([np.percentile(x, 90) for x in fcmatrix])

# Threshold each row of the matrix by setting values below 90th percentile to 0
for i in range(fcmatrix.shape[0]):
  print("Row %d" % i)
  fcmatrix[i, fcmatrix[i,:] < perc[i]] = 0

# Check for minimum value
print("Minimum value is %f" % fcmatrix.min())

# how many nodes have negative values
# Count negative values per row
neg_values = np.array([sum(fcmatrix[i,:] < 0) for i in range(N)])
print("Negative values occur in %d rows" % sum(neg_values > 0))


# example subject 1 had no negative values survive, but imagine other subjects might we set these to zero
fcmatrix[fcmatrix < 0] = 0

# Now we are dealing with sparse vectors. Cosine similarity is used as affinity metric
aff = 1 - pairwise_distances(fcmatrix, metric = 'cosine')


### Commented out saving this: cubic can't hang with multiple 17734x17734 matrices saved per subject (few GB each)
# Save affinity matrix
# savepath= "/cbica/projects/pinesParcels/data/CombinedData/" + str(sid) + "/vertexwise_cos_affinmat.npy"

# np.save(savepath, aff)

# save checkpoint reached, now compute dmap
from mapalign import embed
emb, res = embed.compute_diffusion_map(aff, alpha = 0.5, return_result=True)

savepathe= "/cbica/projects/pinesParcels/data/CombinedData/" + str(sid) + "/vertexwise_emb.npy"
savepathr= "/cbica/projects/pinesParcels/data/CombinedData/" + str(sid) + "/vertexwise_res.npy"

np.save(savepathe, emb)
np.save(savepathr, res)

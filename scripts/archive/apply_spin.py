### apply spins
import nibabel as nb
import scipy
from scipy.io import loadmat
import sys
# first arg is subj
subj = sys.argv[1]
# second arg is spin #
spinNum = int(sys.argv[2])
# load in PG templates
GradsL=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_L_10k.func.gii')
GradsR=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients_R_10k.func.gii')
# extract 1st grad
PGL=GradsL.darrays[0]
PGR=GradsR.darrays[0]
# load in spinmat
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
spinmatFn=parentfp + str(subj) + '_spunions.mat'
spinmat=loadmat(spinmatFn)
# left
Lrot=spinmat['bigrotl']
# right
Rrot=spinmat['bigrotr']
# apply spinNum S (s - 1 because python indexing is dumb)
PGL.data=Lrot[(spinNum-1),:]
PGR.data=Rrot[(spinNum-1),:]
GradsL.darrays[0]=PGL
GradsR.darrays[0]=PGR
# save em
Lfn=parentfp + str(subj) + '_spunPGL.func.gii'
Rfn=parentfp + str(subj) + '_spunPGR.func.gii'
nb.save(GradsL,Lfn)
nb.save(GradsR,Rfn)




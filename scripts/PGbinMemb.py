import nibabel as nb
import hcp_utils as hcp
import numpy as np

# load in principal gradient
PG=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii')

# load in groupCons
GC=nb.load('/gpfs/fs001/cbica/projects/abcdfnets/results/SingleParcellation/RobustInitialization_Cifti_Surf/robust_initVHP.dscalar.nii')

# initialize empty array for gradient bins
MembByBins=np.zeros((2971,20))

# extract cortical only from both
PGc=PG.dataobj[:,hcp.struct.cortex]
GCc=GC.dataobj[:,hcp.struct.cortex]

for b in range(20):
	gradPrctile=np.percentile(PGc[0,:],(b*5))
	gradPrctile_upper=np.percentile(PGc[0,:],(b*5)+5)
	print(b)
	print(gradPrctile)
	print(gradPrctile_upper)
	boolean_of_interest=np.logical_and(PGc[0,:] > gradPrctile, PGc[0,:] < gradPrctile_upper)
	PGindices=np.nonzero(boolean_of_interest)
	PGindices_array=np.array(PGindices)
	#MembByBins[0:(PGindices_array.shape[1]),b]=PGc[0,PGindices_array]
	MembByBins[0:(PGindices_array.shape[1]),b]=GCc[0,PGindices_array]

	

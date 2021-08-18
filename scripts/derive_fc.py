# imports
import scipy
import sys
import nibabel as nb
import hcp_utils as hcp
import numpy as np

# all the scan types
tasks=['rest','SST','nback','mid']

# Subject is set to the passed argument
subj = sys.argv[1]

# load individualized parcellation
parcelLoc='/cbica/projects/abcdfnets/results/SingleParcel_1by1/' + str(subj) + '/' + str(subj) + '_Parcel.dscalar.nii'
parcel=nb.load(parcelLoc)
parcelCort=parcel.dataobj[:,hcp.struct.cortex]

# parent filepath for time series
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'

# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'

# initialize master time series for aggregate fc matrix for dmap
aggregateTS=np.empty([1,59412])

# for each task
for T in range(len(tasks)):
	# load Concat-maskedtime series
	filepath=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '.dtseries.nii'
	subjData=nb.load(filepath)
	# extract cortical TS
	subjDataCort=subjData.dataobj[:,hcp.struct.cortex]
	# stack on to aggregate TS
	aggregateTS=np.append(aggregateTS,subjDataCort,axis=0)
	# gen fc matrix (All)
	GO_fcMat=np.corrcoef(subjDataCort,rowvar=False)	
	# initialize network-level FC matrix
	NL_fcMat=np.zeros((17,17))
	# Individual network metrics (FC matrix of all 17 networks, summary-level)
	for N in range(17):
		# index of which vertices belong to this network
		# +1 to resolve python-everything else counting discordance
		_,NetIndex=np.where(parcelCort==(N+1))
		# index all other networks
		for OtherNet in range(17):
			_,OtherNetIndex=np.where(parcelCort==(OtherNet+1))
			# get elements from both, np.newaxis makes them play nicely together
			NetOtherNetFC=np.mean(GO_fcMat[NetIndex[:,np.newaxis],OtherNetIndex])
			NL_fcMat[N,OtherNet]=NetOtherNetFC
		# correct for same-cell autocorrelation: 
		# subtract ((Vertex#inNet*.5)/NumberOfCells) from mean
		# .5 because diagonals only show up once in symmetric fc, other elements show up tiwce
		WithinCon=NL_fcMat[N,N]
		WithinConAdjusted=WithinCon-((len(NetIndex)*.5)/((len(NetIndex)*len(NetIndex))))
		NL_fcMat[N,N]=WithinConAdjusted
	# save out 17x17 for within-task network stats
	NLfp=childfp+str(subj)+'_NLFcmat'+tasks[T]
	np.save(NLfp,NL_fcMat)

# remove null column of aggregate TS
aggregateTS=np.delete(aggregateTS,0,0)
# make an all-TR fc matrix
aggregateFC=np.corrcoef(aggregateTS,rowvar=False)
# save out full grayord fc mat for diffusion map embedding
Gfp=parentfp+str(subj)+'_FullFCmat'
np.save(Gfp,aggregateFC)

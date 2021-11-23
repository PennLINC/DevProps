# imports
import scipy
import sys
import nibabel as nb
import hcp_utils as hcp
import numpy as np
import os

# function from:
# https://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)
    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

# all the scan types, range(1) only probes rest
tasks=['rest','SST','nback','mid']

# Subject is set to the passed argument
subj = sys.argv[1]

# load individualized parcellation
# parcelLoc='/cbica/projects/abcdfnets/results/SingleParcel_1by1/' + str(subj) + '/' + str(subj) + '_Parcel.dscalar.nii'
# parcel=nb.load(parcelLoc)
# parcelCort=parcel.dataobj[:,hcp.struct.cortex]

# parent filepath for time series
parentfp='/cbica/projects/pinesParcels/results/PWs/PreProc/' + str(subj) + '/'

# child (output) filepath is equivalent to parent in this instance
childfp=parentfp;

# initialize master time series for aggregate fc matrix for dmap
aggregateTS=np.empty([1,59412])

# for each task
# Oct 26 2021: updated to range (1) to ONLY run on resting state: still printed out as AggTS 
for T in range(1):
	# load Concat-maskedtime series
	filepath=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '.dtseries.nii'
	# if file exists
	if os.path.exists(filepath):
		subjData=nb.load(filepath)
		# extract cortical TS
		subjDataCort=subjData.dataobj[:,hcp.struct.cortex]
		# stack on to aggregate TS
		aggregateTS=np.append(aggregateTS,subjDataCort,axis=0)
		# initialize network-level FC matrix
		NL_fcMat=np.zeros((17,17))
		# Individual network metrics (FC matrix of all 17 networks, summary-level)
		for N in range(17):
			# index of which vertices belong to this network
			# +1 to resolve python-everything else counting discordance
			_,NetIndex=np.where(parcelCort==(N+1))
			NetTS=subjDataCort[:,NetIndex]
			# index all other networks
			for OtherNet in range(17):
				# initialize array of Net * OtherNet corrs (vertexwise)
				_,OtherNetIndex=np.where(parcelCort==(OtherNet+1))
				OtherNetTS=subjDataCort[:,OtherNetIndex]
				# correlations between Net and OtherNet elements
				Cors=corr2_coeff(NetTS.T,OtherNetTS.T)
				NL_fcMat[N,OtherNet]=np.mean(Cors)
			# correct for same-cell autocorrelation: 
			# subtract ((Vertex#inNet*.5)/NumberOfCells) from mean
			# .5 because diagonals only show up once in symmetric fc, other elements show up twice
			WithinCon=NL_fcMat[N,N]
			WithinConAdjusted=WithinCon-((len(NetIndex)*.5)/((len(NetIndex)*len(NetIndex))))
			NL_fcMat[N,N]=WithinConAdjusted
		# save out 17x17 for within-task network stats
		NLfp=childfp+str(subj)+'_NLFcmat'+tasks[T]
		np.save(NLfp,NL_fcMat)

# remove null column of aggregate TS
aggregateTS=np.delete(aggregateTS,0,0)
# create new image axis to represent expanded TRs in merged file
newAxis=nb.cifti2.SeriesAxis(start=0,size=aggregateTS.shape[0],step=1)
# create new scalars to populate non-cort with 0, cort with time series
FullBrain=np.zeros((aggregateTS.shape[0],91282))
# make cortical values non-zero (aggregated time series)
FullBrain[:,hcp.struct.cortex]=aggregateTS
# preserved CortgrayOrd axis
gOrdAx=subjData.header.get_axis(1)
new_img=nb.Cifti2Image(FullBrain,(newAxis,gOrdAx))
# save out aggregate time series as cifti for downsampling and diffusion map embedding
subjAggTS=parentfp+str(subj)+'_AggTS.dtseries.nii'
new_img.to_filename(subjAggTS)

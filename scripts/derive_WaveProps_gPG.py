# Normalization, binning of vertexwise time series, wave property saveout 
import scipy
import nibabel as nb
import numpy as np
from scipy import stats
from scipy import signal
from scipy.signal import find_peaks
from scipy.signal import butter
from scipy import signal
from numpy import genfromtxt
import sys
import sklearn
from sklearn import linear_model
import hcp_utils as hcp

#### load in nonSubj TS data
# load in principal gradient
PG=nb.load('/cbica/projects/abcdfnets/data/hcp.gradients.dscalar.nii')

# Subject is set to the passed argument
subj = sys.argv[1]

# all the scan types
tasks=['rest','SST','nback','mid'];

# parent filepath
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# for each task
for T in range(len(tasks)):
	# initialize big array for distribution of all PG delay correlations
	PGD_arr=[]
	# load in time series (masked, bp filtered)
	filepath=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '.dtseries.nii'
	subjData=nb.load(filepath)
	# select just cortex
	subjDataCort=subjData.dataobj[:,hcp.struct.cortex]
	PGCort=PG.dataobj[:,hcp.struct.cortex]
	# load in time series
	procTS=subjDataCort
	# normalize to mean and SD
	Avg=np.mean(procTS,axis=0)
	SD=np.std(procTS,axis=0)
	procTS=(procTS-Avg)/SD
	# convert TS to numpy array to allow for indexing
	procTS=np.array(procTS)
	# initialize empty array for gradient bins
	procTS_bins=np.zeros((len(procTS),70))
	# bin vertices at same position on gradient
	for b in range(70):
		gradPrctile=np.percentile(PGCort[0,:],(b*1.42857))
		gradPrctile_upper=np.percentile(PGCort[0,:],(b*1.42857)+1.42857)
		# index of vertices belonging to this percentile
		boolean_of_interest=np.logical_and(PGCort[0,:] > gradPrctile, PGCort[0,:] < gradPrctile_upper)
		PGindices=np.nonzero(boolean_of_interest)
		PGindices_array=np.array(PGindices)
		# initialize array of all these vertices to average over
		initPGbin=procTS[:,PGindices_array]
		# average signal over this bin 
		meanSig=np.mean(initPGbin,axis=2)
		meanSig=meanSig[:,0]
		# plop into ProcTS_bins
		procTS_bins[:,b]=meanSig
	
	# load in global signal
	GSfFP=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '_GS.csv'
	GS=np.genfromtxt(GSfFP,delimiter=",")
	# normalize GS
	GAvg=np.mean(GS)
	GSD=np.std(GS)
	GS=((GS-GAvg)/GSD)
	# calculate GS troughs with negative find_peaks
	GS_troughs, _ = find_peaks(-GS, distance=8)
	# make an array for each percentile bin and each trough
	troughsNum=len(GS_troughs)-1
	delayMatrix=np.zeros((70,troughsNum))
	# and a relative magnitude matrix
	magMatrix=np.zeros((70,troughsNum))
	# for each trough-trough interval, find peak of bin timeseries
	for t in range(troughsNum):
		tstart=GS_troughs[t]
		tend=GS_troughs[t+1]
		# get GS peak here
		GS_peak, _ = find_peaks(GS[tstart:tend],distance=(tend-tstart))
		for b in range(70):
			# isolate time series sequence
			iso_ts=procTS_bins[tstart:tend,b]
			# find peak in this bin (set min distance to be temporal width of bin)
			peak, _ =find_peaks(iso_ts,distance=(tend-tstart))
			# determine distance from GS peak
			distanceFGSP=peak-GS_peak
			# if peak exists, add to matrix
			if ((len(peak) !=0) and (len(GS_peak) !=0)):
				delayMatrix[b,t]=distanceFGSP
				# record magnitude of normalized signal as point of peak
				magMatrix[b,t]=iso_ts[peak]
			else:
				delayMatrix[b,t]=999
				magMatrix[b,t]=999
	
	# ID columns with < 20% 999s, sep out non-999 values
	# matrix to count instances of no peak detection by PG bin
	npMatrix=np.zeros((70,troughsNum))
	for b in range(70):
		nineninenines = delayMatrix[b,:] == 999
		npMatrix[b,:]=nineninenines
	
	# number of bins w/o detected peak per wave
	noPeakPwave=sum(npMatrix)
	# if peak detected in > 50% of waves, keep it
	mostHavePeaks=delayMatrix[:,noPeakPwave<35]
	# replace 999s with NAs	
	mostHavePeaks[mostHavePeaks==999]=np.nan
	# get nan index for stats
	nas = np.isnan(mostHavePeaks)
	# opposite is valid
	valid = (nas - 1) * -1
	# calculate distribution of correlations of (PG location, trough offset) for each interval
	CorDistr=np.zeros((1,mostHavePeaks.shape[1]))
	# calculate relative magnitude of wave over its course in units of normalized signal
	Wslopes=np.zeros((1,mostHavePeaks.shape[1]))
	# wave speeds (in TRs)
	Wspeeds=np.zeros((1,mostHavePeaks.shape[1]))
	# wave origin point (largest negative delay)
	Worigin=np.zeros((1,mostHavePeaks.shape[1]))
	for i in range(mostHavePeaks.shape[1]):
		# index out nans
		validVec=valid[:,i]!=0
		delayNoNan=mostHavePeaks[validVec,i]
		PGnoNan=np.arange(70)[validVec]
		CorDistr[0,i], _ =stats.spearmanr(delayNoNan,PGnoNan)
		# set x to relative magnitude
		my_x = magMatrix[validVec,i]
		# set y to position in gradient
		slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(my_x, PGnoNan)
		Wslopes[0,i] = slope
		# furthest delayed minus furthest proceeding GS peak
		Wspeeds[0,i] = max(delayNoNan)-min(delayNoNan)
		# where in the gradient is the further proceeding peak?
		val, idx = min((val, idx) for (idx, val) in enumerate(delayNoNan))
		# match to PG with invalids removed
		Worigin[0,i]=PGnoNan[idx]

	# save out GW # x 4 matrix for this subj for this task: PG*delay corr, Speed, slope, and origin for each GW
	saveoutMat=np.zeros((5,mostHavePeaks.shape[1]))
	saveoutMat[0,:]=CorDistr
	saveoutMat[1,:]=Wspeeds
	saveoutMat[2,:]=Worigin
	saveoutMat[3,:]=Wslopes
	# this row will be redudant, will only have one value, number of TRs
	saveoutMat[4,:]=len(GS)
	saveFN=childfp + str(subj) + '_' + str(tasks[T]) + '_waveProps_gPG.csv'
	np.savetxt(saveFN,saveoutMat,delimiter=",")

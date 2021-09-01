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
import os.path
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
# Subject is set to the passed argument
subj = sys.argv[1]
# all the scan types
tasks=['rest','SST','nback','MID'];
# parent filepath
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# load in PG
subjPGfn=childfp+str(subj)+'_PG_LR_32k.dscalar.nii'
PG=nb.load(subjPGfn)
PGdataObject=PG.dataobj
# for each task that exists for this subject
for T in range(len(tasks)):
	# initialize big array for distribution of all PG delay correlations
	PGD_arr=[]
	# load in continuous segment indices
	CSIfp=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-' + tasks[T] + '_ValidSegments.txt'
	# add "if exists" clause to skip processing attempts on non-existent scans
	if os.path.isfile(CSIfp): 
		CSI=np.genfromtxt(CSIfp,delimiter=',')
		# load in time series (masked, bp filtered)
		filepath=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '.dtseries.nii'
		subjData=nb.load(filepath)
		# select just cortex
		subjDataCort=subjData.dataobj[:,hcp.struct.cortex]
		PGCort=PGdataObject[:,hcp.struct.cortex]
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
			# 1.42857 is modifier to allow "70" to fill 70 percentile bins all the way up to 100
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
		# get number of valid segments
		SegShape=CSI.shape
		SegNum=SegShape[0]
		# initialize delay and magnitude matrix (unknown how many wave instances there will be at this point, bin # known
		delayMatrix=np.zeros((70,1))
		magMatrix=np.zeros((70,1))
		# and signal matrix to plot wave unfolding in dif. pg bins - 100 is erring on the side of inclusion
		sigMatrix=np.zeros((100,6))
		totalTroughNum=0
		# for each continuous segment
		for seg in range(SegNum):
			SegStart=CSI[seg,0]
			# python starts at 0, matlab starts at 1
			SegStartInd=int(SegStart-1)
			SegSpan=int(CSI[seg,1])
			# Segment span accounts for first frame, so adding them for indexing is too inclusive by 1
			# -1 to duration for start + duration indexing
			GSinSeg=GS[SegStartInd:(SegStartInd+SegSpan-1)]
			# get time series of grayOrds in this segment
			procTS_bins_inSeg=procTS_bins[SegStartInd:(SegStartInd+SegSpan-1)]
			# calculate GS troughs with negative find_peaks
			GS_troughs, _ = find_peaks(-GSinSeg, distance=8)
			# if there are at least two troughs, we can look at delay in the peak b/w them
			if len(GS_troughs) > 1:
				# make an array for each percentile bin and each trough
				troughsNum=len(GS_troughs)-1
				delayMatrix_Seg=np.zeros((70,troughsNum))
				# and a relative magnitude matrix
				magMatrix_Seg=np.zeros((70,troughsNum))
				# for each trough-trough interval, find peak of bin timeseries
				for t in range(troughsNum):
					tstart=GS_troughs[t]
					tend=GS_troughs[t+1]
					# get GS peak here
					GS_peak, _ = find_peaks(GSinSeg[tstart:tend],distance=(tend-tstart))
					for b in range(70):
						# isolate time series sequence
						iso_ts=procTS_bins_inSeg[tstart:tend,b]
						# find peak in this bin (set min distance to be temporal width of bin)
						peak, _ =find_peaks(iso_ts,distance=(tend-tstart))
						# determine distance from GS peak
						distanceFGSP=peak-GS_peak
						# if peak exists, add to matrix
						if ((len(peak) !=0) and (len(GS_peak) !=0)):
							delayMatrix_Seg[b,t]=distanceFGSP
							# record magnitude of normalized signal as point of peak
							magMatrix_Seg[b,t]=iso_ts[peak]
						else:
							delayMatrix_Seg[b,t]=999
							magMatrix_Seg[b,t]=999
					# record the GS across this wave instance in plot signal matrix
					thisWaveSigMatrix=np.zeros((100,6))
					# record from evenly spaced pgbins as well
					thisWaveSigMatrix[0:(tend-tstart),0]=GSinSeg[tstart:tend]
					thisWaveSigMatrix[0:(tend-tstart),0]=GSinSeg[tstart:tend]
					thisWaveSigMatrix[0:(tend-tstart),1]=procTS_bins_inSeg[tstart:tend,0]
					thisWaveSigMatrix[0:(tend-tstart),2]=procTS_bins_inSeg[tstart:tend,17]
					thisWaveSigMatrix[0:(tend-tstart),3]=procTS_bins_inSeg[tstart:tend,35]
					thisWaveSigMatrix[0:(tend-tstart),4]=procTS_bins_inSeg[tstart:tend,53]
					thisWaveSigMatrix[0:(tend-tstart),5]=procTS_bins_inSeg[tstart:tend,69]
					# tag it onto the master sigMatrix (append into 3d, index out later)
					sigMatrix=np.dstack((sigMatrix,thisWaveSigMatrix))			
				delayMatrix=np.concatenate((delayMatrix,delayMatrix_Seg),axis=1)
				magMatrix=np.concatenate((magMatrix,magMatrix_Seg),axis=1)
				totalTroughNum += troughsNum
		# remove initialization volume of delay and mag matrices
		delayMatrix=delayMatrix[:,1:]
		magMatrix=magMatrix[:,1:]
		sigMatrix=sigMatrix[:,:,1:]
		# ID columns with < 20% 999s, sep out non-999 values
		# matrix to count instances of no peak detection by PG bin
		npMatrix=np.zeros((70,totalTroughNum))
		for b in range(70):
			nineninenines = delayMatrix[b,:] == 999
			npMatrix[b,:]=nineninenines
		# number of bins w/o detected peak per wave
		noPeakPwave=sum(npMatrix)
		# if peak detected in most PG bins, keep it
		mostHavePeaks=delayMatrix[:,noPeakPwave<35]
		# replace 999s with NAs	
		mostHavePeaks[mostHavePeaks==999]=np.nan
		# use same thresholding for sigMatrix
		sigMatrix=sigMatrix[:,:,noPeakPwave<35]
		# for surviving waves
		for m in range(mostHavePeaks.shape[1]):
			plotGS=sigMatrix[:,0,m]
			plt.plot(plotGS[np.nonzero(plotGS)],c='black')
			plotPGB1=sigMatrix[:,1,m]
			plt.plot(plotPGB1[np.nonzero(plotGS)],c='#070291')
			plotPGB2=sigMatrix[:,2,m]
			plt.plot(plotPGB2[np.nonzero(plotGS)],c='#8202ac')
			plotPGB3=sigMatrix[:,3,m]
			plt.plot(plotPGB3[np.nonzero(plotGS)],c='#c8016a')
			plotPGB4=sigMatrix[:,4,m]
			plt.plot(plotPGB4[np.nonzero(plotGS)],c='#e32b01')
			plotPGB5=sigMatrix[:,5,m]
			plt.plot(plotPGB5[np.nonzero(plotGS)],c='#ffe700')
			figName=childfp+str(subj)+'_'+str(tasks[t])+'_Wave'+str(m)+'.png'
			plt.savefig(figName,bbox_inches='tight')
			plt.close()
		# print out wave instances as pyplot
		for m in range(mostHavePeaks.shape[1]):
			plt.plot(mostHavePeaks[:,m]);
			figName=childfp+str(subj)+'_'+str(tasks[t])+'_Delay'+str(m)+'.png'
			plt.savefig(figName,bbox_inches='tight')
			plt.close()
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
		saveFN=childfp + str(subj) + '_' + str(tasks[T]) + '_waveProps.csv'
		np.savetxt(saveFN,saveoutMat,delimiter=",")
		# report difference between all instances of PW and those not meeting >80% threshold
		UnThreshThreshDif=delayMatrix.shape[1]-mostHavePeaks.shape[1]
		print('OG delayMat wave count: ' + str(totalTroughNum))
		print('Waves removed w/ > 50% thresh: ' +str(UnThreshThreshDif))
		saveFN_thr=childfp + str(subj) + '_' + str(tasks[T]) + '_ThreshedWaves'
		np.savetxt(saveFN_thr,[UnThreshThreshDif])

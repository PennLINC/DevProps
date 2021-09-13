# basis time series-based wave property saveout 
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
import mat73
# Subject is set to the passed argument
subj = sys.argv[1]
# all the scan types
tasks=['rest','SST','nback','MID'];
# parent filepath
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# load in Basis time series
TSLoc='/cbica/projects/abcdfnets/results/SingleParcel_1by1/' + str(subj) + '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat'
UV=mat73.loadmat(TSLoc)
# extract TRs by networks matrix
TS=UV['U'][0]
# initialize counter to measure when last task left off and where new begins in merged TS
prevEnd=0
# load in PG
subjPGfn=childfp+str(subj)+'_PG_LR_32k.dscalar.nii'
PG=nb.load(subjPGfn)
PGdataObject=PG.dataobj
# for each task that exists for this subject
for T in range(len(tasks)):
	# initialize "clean" time series with shorter-than-threshold windows removed
	cTS=np.zeros((3,0))
	# initialize big array for distribution of all PG delay correlations
	PGD_arr=[]
	# load in continuous segment indices
	CSIfp=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-' + tasks[T] + '_ValidSegments_Trunc.txt'
	# add "if exists" clause to skip processing attempts on non-existent scans
	if os.path.isfile(CSIfp): 
		CSI=np.genfromtxt(CSIfp,delimiter=',')
		# load in unthreshed segment indices
		SIfp=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-' + tasks[T] + '_ValidSegments_Unthr.txt'
		SI=np.genfromtxt(SIfp,delimiter=',')
		# get span of this task relative to fully concatenated TS
		startOfThisTask=prevEnd
		# extract length of "passing" TRs from unthr file, -1 at end bc length of window includes startframe
		passingTRs=SI[-1,0]+SI[-1,1]-1
		EndOfThisTask=prevEnd+passingTRs
		# note: end of task is printed in matlabic terms: this is bc python excludes last element in x:y ranges
		print(EndOfThisTask)
		# Extract the task time series
		TS_task=TS[int(startOfThisTask):int(EndOfThisTask),:]
		# now that prevEnd is used, set it for the next iteration
		prevEnd=EndOfThisTask
		# Isolate cont. windows from FullValidSeg files
		fSIfp=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-' + tasks[T] + '_ValidSegments_Full.txt'
		fSI=np.genfromtxt(fSIfp,delimiter=',')
		# match bandpassing script for extraction, tag on valid segments to the initialized "clean" time series
		# matching starts with getting number of segments
		numSegs=fSI.shape[0]
		for s in range(numSegs):
			### always a battle beteween pythonic and matlabic indexing
			# START at index -1, because python starts at 0
			pythonicStart=int(fSI[s,0]-1)
			# BUT, end at same point because python range (x:y) is exclusive
			pythonicEnd=int(pythonicStart+fSI[s,1])
			segment=TS_task[pythonicStart:pythonicEnd,:]
			# Intermed. version. Isolating Aud, Mot, Vis, DAN, VAN, FPN, DMNs
			segAud=segment[:,15]
			segMot=np.mean(segment[:,(12,10,3,1)],axis=1)
			segVis=np.mean(segment[:,(5,9)],axis=1)
			segDAN=np.mean(segment[:,(13,4)],axis=1)
			segVAN=np.mean(segment[:,(8,6)],axis=1)
			segFPN=np.mean(segment[:,(2,14,16)],axis=1)
			segDMN=np.mean(segment[:,(0,7,11)],axis=1)
			stacked=np.vstack((segAud,segMot,segVis,segDAN,segVAN,segFPN,segDMN))
			cTS=np.concatenate((cTS,stacked),axis=1)
		# transpose because rest of this script follows that formatting
		cTS=np.transpose(cTS)
		# check that procTS matches truncated valid segments file
		numTRsPTS=cTS.shape[0]
		numTRsVS=CSI[-1,0]+CSI[-1,1]-1
		if numTRsPTS != numTRsVS:
			print(numTRsPTS)
			print(numTRsVS)
			raise Exception('TRs from Valid Segments txt and cifti do not match')
		# convert TS to numpy array to allow for indexing
		procTS=np.array(cTS)
		# load in global signal
		GSfFP=parentfp + str(subj) + '_p2mm_masked_filtered_' + tasks[T] + '_GS.csv'
		GS=np.genfromtxt(GSfFP,delimiter=",")
		if GS.shape[0] != numTRsVS:
			print(GS.shape[0])
			print(numTRsVS)
			raise Exception('TRs from Valid Segments txt and GS do not match')
		# normalize GS
		GAvg=np.mean(GS)
		GSD=np.std(GS)
		GS=((GS-GAvg)/GSD)
		# get number of valid segments
		SegShape=CSI.shape
		SegNum=SegShape[0]
		# initialize delay and magnitude matrix (unknown how many wave instances there will be at this point, bin # known
		delayMatrix=np.zeros((3,1))
		magMatrix=np.zeros((3,1))
		# and signal matrix to plot wave unfolding in dif. pg bins - 100 is erring on the side of inclusion
		sigMatrix=np.zeros((100,3))
		totalTroughNum=0
		# and an empty 2-column matrix to record which segment and how far into segment detected GS troughs occur
		waveTRs=np.zeros((1,3))
		# initialize peak vector
		peakVec=[]
		# initialize TR tracker number
		trTrackVec=[]
		# for each continuous segment
		for seg in range(SegNum):
			SegStart=CSI[seg,0]
			# python starts at 0, matlab starts at 1
			SegStartInd=int(SegStart-1)
			SegSpan=int(CSI[seg,1])
			# Segment span accounts for first frame, so adding them for indexing is too inclusive by 1
			# -1 to duration for start + duration indexing
			GSinSeg=GS[SegStartInd:(SegStartInd+SegSpan-1)]
			# get basis time series in this segment
			procTS_bins_inSeg=procTS[SegStartInd:(SegStartInd+SegSpan-1),:]
			# find peaks in Aud,Mot,Vis,DAN,VAN,FPN,DMN
			Gpeaks,_=find_peaks(GSinSeg,distance=8)
			Audpeaks,_=find_peaks(procTS_bins_inSeg[:,0],distance=8)
			Motpeaks,_=find_peaks(procTS_bins_inSeg[:,1],distance=8)
			Vispeaks,_=find_peaks(procTS_bins_inSeg[:,2],distance=8)
			DANpeaks,_=find_peaks(procTS_bins_inSeg[:,3],distance=8)
			VANpeaks,_=find_peaks(procTS_bins_inSeg[:,4],distance=8)
			FPNpeaks,_=find_peaks(procTS_bins_inSeg[:,5],distance=8)
			DMNpeaks,_=find_peaks(procTS_bins_inSeg[:,6],distance=8)
			# for each TR
			for TR in range(SegSpan):
				# is there a GS peak?
				if TR in Gpeaks:
					peakVec=np.append(peakVec,0)
					trTrackVec=np.append(trTrackVec,TR)
				# is there an Aud peak?
				if TR in Audpeaks:
					peakVec=np.append(peakVec,1)
					trTrackVec=np.append(trTrackVec,TR)
				# is there a Mot peak?
				if TR in Motpeaks:
					peakVec=np.append(peakVec,2)
					trTrackVec=np.append(trTrackVec,TR)
				# is there a Vis peak?
				if TR in Vispeaks:
					peakVec=np.append(peakVec,3)
					trTrackVec=np.append(trTrackVec,TR)
				if TR in DANpeaks:
					peakVec=np.append(peakVec,4)
					trTrackVec=np.append(trTrackVec,TR)
				if TR in VANpeaks:
					peakVec=np.append(peakVec,5)
					trTrackVec=np.append(trTrackVec,TR)
				if TR in FPNpeaks:
					peakVec=np.append(peakVec,6)
					trTrackVec=np.append(trTrackVec,TR)
				if TR in DMNpeaks:
					peakVec=np.append(peakVec,7)
					trTrackVec=np.append(trTrackVec,TR)
			# append both with end of segment marker
			peakVec=np.append(peakVec,-1)
			trTrackVec=np.append(trTrackVec,-1)
		# save out peak vec and TRvec
		PsaveFN=childfp + str(subj) + '_' + str(tasks[T]) + '_PeakSeq_int'
		np.savetxt(PsaveFN,peakVec)
		TsaveFN=childfp + str(subj) + '_' + str(tasks[T]) + '_TRSeq_int'
		np.savetxt(TsaveFN,trTrackVec)

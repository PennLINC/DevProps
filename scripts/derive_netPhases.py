import nibabel as nb
import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy import signal
import copy
import scipy
import sys
subj = sys.argv[1]

##### load in proc. data
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
filepath=parentfp + str(subj) + '_p2mm_masked_filtered_rest.dtseries.nii'
subjData=nb.load(filepath)
# extract cortex
subjDataCort=subjData.dataobj[:,hcp.struct.cortex]

##### load in parcellation
parcelLoc='/cbica/projects/abcdfnets/results/SingleParcel_1by1/' + str(subj) + '/' + str(subj) + '_Parcel.dscalar.nii'
parcel=nb.load(parcelLoc)
ParcelC=parcel.dataobj[:,hcp.struct.cortex]
# get indices of parcels
DM_inds12=np.where(ParcelC==12)
DM_inds1=np.where(ParcelC==1)
DM_inds8=np.where(ParcelC==8)
DM_inds=np.hstack((DM_inds12,DM_inds1,DM_inds8))
Mot_inds2=np.where(ParcelC==2)
Mot_inds4=np.where(ParcelC==4)
Mot_inds11=np.where(ParcelC==11)
Mot_inds13=np.where(ParcelC==13)
Vis_inds6=np.where(ParcelC==6)
Vis_inds10=np.where(ParcelC==10)
DAN_inds14=np.where(ParcelC==14)
DAN_inds5=np.where(ParcelC==5)
VAN_inds9=np.where(ParcelC==9)
VAN_inds7=np.where(ParcelC==7)
FPN_inds3=np.where(ParcelC==3)
FPN_inds15=np.where(ParcelC==15)
FPN_inds17=np.where(ParcelC==17)
Aud_inds=np.where(ParcelC==16)

######## load in valid segments file
CSIfp=parentfp + str(subj) + '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'
CSI=np.genfromtxt(CSIfp,delimiter=',')
SegShape=CSI.shape
SegNum=SegShape[0]

# initialize cross-segment amp and phase to stack within-seg amp and phase onto
DM1phases_Agg=np.zeros(0)
DM8phases_Agg=np.zeros(0)
DM12phases_Agg=np.zeros(0)
FP3phases_Agg=np.zeros(0)
FP15phases_Agg=np.zeros(0)
FP17phases_Agg=np.zeros(0)
VAN7phases_Agg=np.zeros(0)
VAN9phases_Agg=np.zeros(0)
DAN5phases_Agg=np.zeros(0)
DAN14phases_Agg=np.zeros(0)
Mot2phases_Agg=np.zeros(0)
Mot4phases_Agg=np.zeros(0)
Mot11phases_Agg=np.zeros(0)
Mot13phases_Agg=np.zeros(0)
Vis6phases_Agg=np.zeros(0)
Vis10phases_Agg=np.zeros(0)
Aud16phases_Agg=np.zeros(0)

DM1Amp_Agg=np.zeros(0)
DM8Amp_Agg=np.zeros(0)
DM12Amp_Agg=np.zeros(0)
FP3Amp_Agg=np.zeros(0)
FP15Amp_Agg=np.zeros(0)
FP17Amp_Agg=np.zeros(0)
VAN7Amp_Agg=np.zeros(0)
VAN9Amp_Agg=np.zeros(0)
DAN5Amp_Agg=np.zeros(0)
DAN14Amp_Agg=np.zeros(0)
Mot2Amp_Agg=np.zeros(0)
Mot4Amp_Agg=np.zeros(0)
Mot11Amp_Agg=np.zeros(0)
Mot13Amp_Agg=np.zeros(0)
Vis6Amp_Agg=np.zeros(0)
Vis10Amp_Agg=np.zeros(0)
Aud16Amp_Agg=np.zeros(0)

# same for DM peak span
PeakIntervalSpans=np.zeros(0)

######## for each segment
for seg in range(SegNum):
	SegStart=CSI[seg,0]
	# python starts at 0, matlab starts at 1
	SegStartInd=int(SegStart-1)
	SegSpan=int(CSI[seg,1])
	# Segment span accounts for first frame, so adding them for indexing is too inclusive by 1
	# -1 to duration for start + duration indexing
	subjDataCor_inSeg=subjDataCor[SegStartInd:(SegStartInd+SegSpan-1),:]
	# get time series assoc. with each parcel
	# MEAN DM
	meanDM=np.mean(subjDataCor_inSeg[:,DM_inds],axis=2)
	DM=meanDM[:,0]
	# DM
	meanDM1=np.mean(subjDataCor_inSeg[:,DM_inds1],axis=2)
	DM1=meanDM1[:,0]
	meanDM8=np.mean(subjDataCor_inSeg[:,DM_inds8],axis=2)
	DM8=meanDM8[:,0]
	meanDM12=np.mean(subjDataCor_inSeg[:,DM_inds12],axis=2)
	DM12=meanDM12[:,0]
	# FP
	meanFPN3=np.mean(subjDataCor_inSeg[:,FPN_inds3],axis=2)
	FP3=meanFPN3[:,0]
	meanFPN15=np.mean(subjDataCor_inSeg[:,FPN_inds15],axis=2)
	FP15=meanFPN15[:,0]
	meanFPN17=np.mean(subjDataCor_inSeg[:,FPN_inds17],axis=2)
	FP17=meanFPN17[:,0]
	# VAN
	meanVAN7=np.mean(subjDataCor_inSeg[:,VAN_inds7],axis=2)
	VAN7=meanVAN7[:,0]
	meanVAN9=np.mean(subjDataCor_inSeg[:,VAN_inds9],axis=2)
	VAN9=meanVAN9[:,0]
	# DAN
	meanDAN5=np.mean(subjDataCor_inSeg[:,DAN_inds5],axis=2)
	DAN5=meanDAN5[:,0]
	meanDAN14=np.mean(subjDataCor_inSeg[:,DAN_inds14],axis=2)
	DAN14=meanDAN14[:,0]
	# MOT
	meanMot2=np.mean(subjDataCor_inSeg[:,Mot_inds2],axis=2)
	Mot2=meanMot2[:,0]
	meanMot4=np.mean(subjDataCor_inSeg[:,Mot_inds4],axis=2)
	Mot4=meanMot4[:,0]
	meanMot11=np.mean(subjDataCor_inSeg[:,Mot_inds11],axis=2)
	Mot11=meanMot11[:,0]
	meanMot13=np.mean(subjDataCor_inSeg[:,Mot_inds13],axis=2)
	Mot13=meanMot13[:,0]
	# VIS
	meanVis6=np.mean(subjDataCor_inSeg[:,Vis_inds6],axis=2)
	Vis6=meanVis6[:,0]
	meanVis10=np.mean(subjDataCor_inSeg[:,Vis_inds10],axis=2)
	Vis10=meanVis10[:,0]
	# AUD
	meanAud=np.mean(subjDataCor_inSeg[:,Aud_inds],axis=2)
	AUD16=meanAud[:,0]
	# get peaks in mean DM
	negpeaks,_=scipy.signal.find_peaks(-DM,distance=10)
	peaks,_=scipy.signal.find_peaks(DM,distance=10)
	# gating logic for peaks
	peaks=peaks[DM[peaks]>0]
	# use to extract interpeak intervals
	intervals=np.diff(peaks)
	IntervalBoolean=np.ones(len(intervals))
	# remove last peak for combination
	cyclicPeaks=peaks[:-1]
	# combine
	PeaksandDurs=np.c_[cyclicPeaks,intervals]
	##### initialize phase, mag, and peak-peak span vectors
	DM1phases=np.zeros(len(PeaksandDurs))
	DM8phases=np.zeros(len(PeaksandDurs))
	DM12phases=np.zeros(len(PeaksandDurs))
	FP3phases=np.zeros(len(PeaksandDurs))
	FP15phases=np.zeros(len(PeaksandDurs))
	FP17phases=np.zeros(len(PeaksandDurs))
	VAN7phases=np.zeros(len(PeaksandDurs))
	VAN9phases=np.zeros(len(PeaksandDurs))
	DAN5phases=np.zeros(len(PeaksandDurs))
	DAN14phases=np.zeros(len(PeaksandDurs))
	Mot2phases=np.zeros(len(PeaksandDurs))
	Mot4phases=np.zeros(len(PeaksandDurs))
	Mot11phases=np.zeros(len(PeaksandDurs))
	Mot13phases=np.zeros(len(PeaksandDurs))
	Vis6phases=np.zeros(len(PeaksandDurs))
	Vis10phases=np.zeros(len(PeaksandDurs))
	Aud16phases=np.zeros(len(PeaksandDurs))

	DM1Amp=np.zeros(len(PeaksandDurs))
	DM8Amp=np.zeros(len(PeaksandDurs))
	DM12Amp=np.zeros(len(PeaksandDurs))
	FP3Amp=np.zeros(len(PeaksandDurs))
	FP15Amp=np.zeros(len(PeaksandDurs))
	FP17Amp=np.zeros(len(PeaksandDurs))
	VAN7Amp=np.zeros(len(PeaksandDurs))
	VAN9Amp=np.zeros(len(PeaksandDurs))
	DAN5Amp=np.zeros(len(PeaksandDurs))
	DAN14Amp=np.zeros(len(PeaksandDurs))
	Mot2Amp=np.zeros(len(PeaksandDurs))
	Mot4Amp=np.zeros(len(PeaksandDurs))
	Mot11Amp=np.zeros(len(PeaksandDurs))
	Mot13Amp=np.zeros(len(PeaksandDurs))
	Vis6Amp=np.zeros(len(PeaksandDurs))
	Vis10Amp=np.zeros(len(PeaksandDurs))
	Aud16Amp=np.zeros(len(PeaksandDurs))

	##### extract network-specific peaks w/r/t average DMN peaks
	for p in range(len(PeaksandDurs)):
		IntervalSpan=PeaksandDurs[p,1]
		# need to pad on both sides to allow for possibility that peak is concurrent w/ dmn peak
		DM1InInterval=DM1[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		DM8InInterval=DM8[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		DM12InInterval=DM12[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		FP3InInterval=FP3[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		FP15InInterval=FP15[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		FP17InInterval=FP17[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		VAN7InInterval=VAN7[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		VAN9InInterval=VAN9[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		DAN5InInterval=DAN5[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		DAN14InInterval=DAN14[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Mot2InInterval=Mot2[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Mot4InInterval=Mot4[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Mot11InInterval=Mot11[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Mot13InInterval=Mot13[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Vis6InInterval=Vis6[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Vis10InInterval=Vis10[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		Aud16InInterval=AUD16[(PeaksandDurs[p,0]-1):(PeaksandDurs[p,0]+PeaksandDurs[p,1]+1)]
		# mark peaks in DM phase, assuming startpoint is 0 and endpoint is 360
		DM1Peak=signal.find_peaks(DM1InInterval,distance=IntervalSpan)[0]
		DM8Peak=signal.find_peaks(DM8InInterval,distance=IntervalSpan)[0]
		DM12Peak=signal.find_peaks(DM12InInterval,distance=IntervalSpan)[0]
		FP3Peak=signal.find_peaks(FP3InInterval,distance=IntervalSpan)[0]
		FP15Peak=signal.find_peaks(FP15InInterval,distance=IntervalSpan)[0]
		FP17Peak=signal.find_peaks(FP17InInterval,distance=IntervalSpan)[0]
		VAN7Peak=signal.find_peaks(VAN7InInterval,distance=IntervalSpan)[0]
		VAN9Peak=signal.find_peaks(VAN9InInterval,distance=IntervalSpan)[0]
		DAN5Peak=signal.find_peaks(DAN5InInterval,distance=IntervalSpan)[0]
		DAN14Peak=signal.find_peaks(DAN14InInterval,distance=IntervalSpan)[0]
		mot2Peak=signal.find_peaks(Mot2InInterval,distance=IntervalSpan)[0]
		mot4Peak=signal.find_peaks(Mot4InInterval,distance=IntervalSpan)[0]
		mot11Peak=signal.find_peaks(Mot11InInterval,distance=IntervalSpan)[0]
		mot13Peak=signal.find_peaks(Mot13InInterval,distance=IntervalSpan)[0]
		vis6Peak=signal.find_peaks(Vis6InInterval,distance=IntervalSpan)[0]
		vis10Peak=signal.find_peaks(Vis10InInterval,distance=IntervalSpan)[0]
		aud16Peak=signal.find_peaks(Aud16InInterval,distance=IntervalSpan)[0]
		# convert to percent through the span (if peaks exist)
		if len(DM1Peak)>0:
		    DM1PeakPerc=DM1Peak/IntervalSpan
		    DM1PeakPercRad=(np.pi*2*DM1PeakPerc)
		    DM1phases[p]=DM1PeakPercRad
		    DM1Amp[p]=DM1InInterval[DM1Peak]
		elif len(DM1Peak)==0:
		    DM1phases[p]=99
		if len(DM8Peak)>0:
		    DM8PeakPerc=DM8Peak/IntervalSpan
		    DM8PeakPercRad=(np.pi*2*DM8PeakPerc)
		    DM8phases[p]=DM8PeakPercRad
		    DM8Amp[p]=DM8InInterval[DM8Peak]
		elif len(DM8Peak)==0:
		    DM8phases[p]=99
		if len(DM12Peak)>0:
		    DM12PeakPerc=DM12Peak/IntervalSpan
		    DM12PeakPercRad=(np.pi*2*DM12PeakPerc)
		    DM12phases[p]=DM12PeakPercRad
		    DM12Amp[p]=DM12InInterval[DM12Peak]
		elif len(DM12Peak)==0:
		    DM12phases[p]=99
		if len(FP3Peak)>0:
		    FP3PeakPerc=FP3Peak/IntervalSpan
		    FP3PeakPercRad=(np.pi*2*FP3PeakPerc)
		    FP3phases[p]=FP3PeakPercRad
		    FP3Amp[p]=FP3InInterval[FP3Peak]
		elif len(FP3Peak)==0:
		    FP3phases[p]=99
		if len(FP15Peak)>0:
		    FP15PeakPerc=FP15Peak/IntervalSpan
		    FP15PeakPercRad=(np.pi*2*FP15PeakPerc)
		    FP15phases[p]=FP15PeakPercRad
		    FP15Amp[p]=FP15InInterval[FP15Peak]
		elif len(FP15Peak)==0:
		    FP15phases[p]=99
		if len(FP17Peak)>0:
		    FP17PeakPerc=FP17Peak/IntervalSpan
		    FP17PeakPercRad=(np.pi*2*FP17PeakPerc)
		    FP17phases[p]=FP17PeakPercRad
		    FP17Amp[p]=FP17InInterval[FP17Peak]
		elif len(FP17Peak)==0:
		    FP17phases[p]=99
		if len(VAN7Peak)>0:
		    VAN7PeakPerc=VAN7Peak/IntervalSpan
		    VAN7PeakPercRad=(np.pi*2*VAN7PeakPerc)
		    VAN7phases[p]=VAN7PeakPercRad
		    VAN7Amp[p]=VAN7InInterval[VAN7Peak]
		elif len(VAN7Peak)==0:
		    VAN7phases[p]=99
		if len(VAN9Peak)>0:
		    VAN9PeakPerc=VAN9Peak/IntervalSpan
		    VAN9PeakPercRad=(np.pi*2*VAN9PeakPerc)
		    VAN9phases[p]=VAN9PeakPercRad
		    VAN9Amp[p]=VAN9InInterval[VAN9Peak]
		elif len(VAN9Peak)==0:
		    VAN9phases[p]=99
		if len(DAN5Peak)>0:
		    DAN5PeakPerc=DAN5Peak/IntervalSpan
		    DAN5PeakPercRad=(np.pi*2*DAN5PeakPerc)
		    DAN5phases[p]=DAN5PeakPercRad
		    DAN5Amp[p]=DAN5InInterval[DAN5Peak]
		elif len(DAN5Peak)==0:
		    DAN5phases[p]=99
		if len(DAN14Peak)>0:
		    DAN14PeakPerc=DAN14Peak/IntervalSpan
		    DAN14PeakPercRad=(np.pi*2*DAN14PeakPerc)
		    DAN14phases[p]=DAN14PeakPercRad
		    DAN14Amp[p]=DAN14InInterval[DAN14Peak]
		elif len(DAN14Peak)==0:
		    DAN14phases[p]=99
		if len(mot2Peak)>0:
		    mot2PeakPerc=mot2Peak/IntervalSpan
		    mot2PeakPercRad=(np.pi*2*mot2PeakPerc)
		    Mot2phases[p]=mot2PeakPercRad
		    Mot2Amp[p]=Mot2InInterval[mot2Peak]
		elif len(mot2Peak)==0:
		    Mot2phases[p]=99
		if len(mot4Peak)>0:
		    mot4PeakPerc=mot4Peak/IntervalSpan
		    mot4PeakPercRad=(np.pi*2*mot4PeakPerc)
		    Mot4phases[p]=mot4PeakPercRad
		    Mot4Amp[p]=Mot4InInterval[mot4Peak]
		elif len(mot4Peak)==0:
		    Mot4phases[p]=99
		if len(mot11Peak)>0:
		    mot11PeakPerc=mot11Peak/IntervalSpan
		    mot11PeakPercRad=(np.pi*2*mot11PeakPerc)
		    Mot11phases[p]=mot11PeakPercRad
		    Mot11Amp[p]=Mot11InInterval[mot11Peak]
		elif len(mot11Peak)==0:
		    Mot11phases[p]=99
		if len(mot13Peak)>0:
		    mot13PeakPerc=mot13Peak/IntervalSpan
		    mot13PeakPercRad=(np.pi*2*mot13PeakPerc)
		    Mot13phases[p]=mot13PeakPercRad
		    Mot13Amp[p]=Mot13InInterval[mot13Peak]
		elif len(mot13Peak)==0:
		    Mot13phases[p]=99
		if len(vis6Peak)>0:
		    vis6PeakPerc=vis6Peak/IntervalSpan
		    vis6PeakPercRad=(np.pi*2*vis6PeakPerc)
		    Vis6phases[p]=vis6PeakPercRad
		    Vis6Amp[p]=Vis6InInterval[vis6Peak]
		elif len(vis6Peak)==0:
		    Vis6phases[p]=99
		if len(vis10Peak)>0:
		    vis10PeakPerc=vis10Peak/IntervalSpan
		    vis10PeakPercRad=(np.pi*2*vis10PeakPerc)
		    Vis10phases[p]=vis10PeakPercRad
		    Vis10Amp[p]=Vis10InInterval[vis10Peak]
		elif len(vis10Peak)==0:
		    Vis10phases[p]=99
		if len(aud16Peak)>0:
		    aud16PeakPerc=aud16Peak/IntervalSpan
		    aud16PeakPercRad=(np.pi*2*aud16PeakPerc)
		    Aud16phases[p]=aud16PeakPercRad
		    Aud16Amp[p]=Aud16InInterval[aud16Peak]
		elif len(aud16Peak)==0:
		    Aud16phases[p]=99

	# add this segment's phases, amps, and spans to the aggregate vectors
	DM1phases_Agg=np.concatenate((DM1phases_Agg, DM1phases), axis=1)
	DM8phases_Agg=np.concatenate((DM8phases_Agg, DM8phases), axis=1)
	DM12phases_Agg=np.concatenate((DM12phases_Agg, DM12phases), axis=1)
	FP3phases_Agg=np.concatenate((FP3phases_Agg, FP3phases), axis=1)
	FP15phases_Agg=np.concatenate((FP15phases_Agg, FP15phases), axis=1)
	FP17phases_Agg=np.concatenate((FP17phases_Agg, FP17phases), axis=1)
	VAN7phases_Agg=np.concatenate((VAN7phases_Agg, VAN7phases), axis=1)
	VAN9phases_Agg=np.concatenate((VAN9phases_Agg, VAN9phases), axis=1)
	DAN5phases_Agg=np.concatenate((DAN5phases_Agg, DAN5phases), axis=1)
	DAN14phases_Agg=np.concatenate((DAN14phases_Agg, DAN14phases), axis=1)
	Mot2phases_Agg=np.concatenate((Mot2phases_Agg, Mot2phases), axis=1)
	Mot4phases_Agg=np.concatenate((Mot4phases_Agg, Mot4phases), axis=1)
	Mot11phases_Agg=np.concatenate((Mot11phases_Agg, Mot11phases), axis=1)
	Mot13phases_Agg=np.concatenate((Mot13phases_Agg, Mot13phases), axis=1)
	Vis6phases_Agg=np.concatenate((Vis6phases_Agg, Vis6phases), axis=1)
	Vis10phases_Agg=np.concatenate((Vis10phases_Agg, Vis10phases), axis=1)
	Aud16phases_Agg=np.concatenate((Aud16phases_Agg, Aud16phases), axis=1)
	# amplitude
	DM1Amp_Agg=np.concatenate((DM1Amp_Agg, DM1Amp), axis=1)
	DM8Amp_Agg=np.concatenate((DM8Amp_Agg, DM8Amp), axis=1)
	DM12Amp_Agg=np.concatenate((DM12Amp_Agg, DM12Amp), axis=1)
	FP3Amp_Agg=np.concatenate((FP3Amp_Agg, FP3Amp), axis=1)
	FP15Amp_Agg=np.concatenate((FP15Amp_Agg, FP15Amp), axis=1)
	FP17Amp_Agg=np.concatenate((FP17Amp_Agg, FP17Amp), axis=1)
	VAN7Amp_Agg=np.concatenate((VAN7Amp_Agg, VAN7Amp), axis=1)
	VAN9Amp_Agg=np.concatenate((VAN9Amp_Agg, VAN9Amp), axis=1)
	DAN5Amp_Agg=np.concatenate((DAN5Amp_Agg, DAN5Amp), axis=1)
	DAN14Amp_Agg=np.concatenate((DAN14Amp_Agg, DAN14Amp), axis=1)
	Mot2Amp_Agg=np.concatenate((Mot2Amp_Agg, Mot2Amp), axis=1)
	Mot4Amp_Agg=np.concatenate((Mot4Amp_Agg, Mot4Amp), axis=1)
	Mot11Amp_Agg=np.concatenate((Mot11Amp_Agg, Mot11Amp), axis=1)
	Mot13Amp_Agg=np.concatenate((Mot13Amp_Agg, Mot13Amp), axis=1)
	Vis6Amp_Agg=np.concatenate((Vis6Amp_Agg, Vis6Amp), axis=1)
	Vis10Amp_Agg=np.concatenate((Vis10Amp_Agg, Vis10Amp), axis=1)
	Aud16Amp_Agg=np.concatenate((Aud16Amp_Agg, Aud16Amp), axis=1)
	# same for DM peak span
	PeakIntervalSpans=np.concatenate((PeakIntervalSpans,PeaksandDurs[:,1]), axis=1)
	
######## end for each segment
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
##### save out phase and magnitude info
DM1pfp=childfp + str(subj) + '_rest_DM1P.csv'
DM1phases_Agg=np.savetxt(DM1phases_Agg, DM1pfp,delimiter=“,”)
DM8pfp=childfp + str(subj) + '_rest_DM8P.csv'
DM8phases_Agg=np.savetxt(DM8phases_Agg, DM8pfp,delimiter=“,”)
DM12pfp=childfp + str(subj) + '_rest_DM12P.csv'
DM12phases_Agg=np.savetxt(DM12phases_Agg, DM12pfp,delimiter=“,”)
FP3pfp=childfp + str(subj) + '_rest_FP3P.csv'
FP3phases_Agg=np.savetxt(FP3phases_Agg, FP3pfp,delimiter=“,”)
FP15pfp=childfp + str(subj) + '_rest_FP15P.csv'
FP15phases_Agg=np.savetxt(FP15phases_Agg, FP15pfp,delimiter=“,”)
FP17pfp=childfp + str(subj) + '_rest_FP17P.csv'
FP17phases_Agg=np.savetxt(FP17phases_Agg, FP17pfp,delimiter=“,”)
VAN7pfp=childfp + str(subj) + '_rest_VAN7P.csv'
VAN7phases_Agg=np.savetxt(VAN7phases_Agg, VAN7pfp,delimiter=“,”)
VAN9pfp=childfp + str(subj) + '_rest_VAN9P.csv'
VAN9phases_Agg=np.savetxt(VAN9phases_Agg, VAN9pfp,delimiter=“,”)
DAN5pfp=childfp + str(subj) + '_rest_DAN5P.csv'
DAN5phases_Agg=np.savetxt(DAN5phases_Agg, DAN5pfp,delimiter=“,”)
DAN14pfp=childfp + str(subj) + '_rest_DAN14P.csv'
DAN14phases_Agg=np.savetxt(DAN14phases_Agg, DAN14pfp,delimiter=“,”)
Mot2pfp=childfp + str(subj) + '_rest_MOT2P.csv'
Mot2phases_Agg=np.savetxt(Mot2phases_Agg, Mot2pfp,delimiter=“,”)
Mot4pfp=childfp + str(subj) + '_rest_MOT4P.csv'
Mot4phases_Agg=np.savetxt(Mot4phases_Agg, Mot4pfp,delimiter=“,”)
Mot11pfp=childfp + str(subj) + '_rest_MOT11P.csv'
Mot11phases_Agg=np.savetxt(Mot11phases_Agg, Mot11pfp,delimiter=“,”)
Mot13pfp=childfp + str(subj) + '_rest_MOT13P.csv'
Mot13phases_Agg=np.savetxt(Mot13phases_Agg, Mot13pfp,delimiter=“,”)
Vis6pfp=childfp + str(subj) + '_rest_VIS6P.csv'
Vis6phases_Agg=np.savetxt(Vis6phases_Agg, Vis6pfp,delimiter=“,”)
Vis10pfp=childfp + str(subj) + '_rest_VIS10P.csv'
Vis10phases_Agg=np.savetxt(Vis10phases_Agg, Vis10pfp,delimiter=“,”)
Aud16pfp=childfp + str(subj) + '_rest_AUD16P.csv'
Aud16phases_Agg=np.savetxt(Aud16phases_Agg, Aud16pfp,delimiter=“,”)

# amplitudes
DM1afp=childfp + str(subj) + '_rest_DM1A.csv'
DM1Amp_Agg=np.savetxt(DM1Amp_Agg, DM1afp,delimiter=“,”)
DM8afp=childfp + str(subj) + '_rest_DM8A.csv'
DM8Amp_Agg=np.savetxt(DM8Amp_Agg, DM8afp,delimiter=“,”)
DM12afp=childfp + str(subj) + '_rest_DM12A.csv'
DM12Amp_Agg=np.savetxt(DM12Amp_Agg, DM12afp,delimiter=“,”)
FP3afp=childfp + str(subj) + '_rest_FP3A.csv'
FP3Amp_Agg=np.savetxt(FP3Amp_Agg, FP3afp,delimiter=“,”)
FP15afp=childfp + str(subj) + '_rest_FP15A.csv'
FP15Amp_Agg=np.savetxt(FP15Amp_Agg, FP15afp,delimiter=“,”)
FP17afp=childfp + str(subj) + '_rest_FP17A.csv'
FP17Amp_Agg=np.savetxt(FP17Amp_Agg, FP17afp,delimiter=“,”)
VAN7afp=childfp + str(subj) + '_rest_VAN7A.csv'
VAN7Amp_Agg=np.savetxt(VAN7Amp_Agg, VAN7afp,delimiter=“,”)
VAN9afp=childfp + str(subj) + '_rest_VAN9A.csv'
VAN9Amp_Agg=np.savetxt(VAN9Amp_Agg, VAN9afp,delimiter=“,”)
DAN5afp=childfp + str(subj) + '_rest_DAN5A.csv'
DAN5Amp_Agg=np.savetxt(DAN5Amp_Agg, DAN5afp,delimiter=“,”)
DAN14afp=childfp + str(subj) + '_rest_DAN14A.csv'
DAN14Amp_Agg=np.savetxt(DAN14Amp_Agg, DAN14afp,delimiter=“,”)
Mot2afp=childfp + str(subj) + '_rest_MOT2A.csv'
Mot2Amp_Agg=np.savetxt(Mot2Amp_Agg, Mot2afp,delimiter=“,”)
Mot4afp=childfp + str(subj) + '_rest_MOT4A.csv'
Mot4Amp_Agg=np.savetxt(Mot4Amp_Agg, Mot4afp,delimiter=“,”)
Mot11afp=childfp + str(subj) + '_rest_MOT11A.csv'
Mot11Amp_Agg=np.savetxt(Mot11Amp_Agg, Mot11afp,delimiter=“,”)
Mot13afp=childfp + str(subj) + '_rest_MOT13A.csv'
Mot13Amp_Agg=np.savetxt(Mot13Amp_Agg, Mot13afp,delimiter=“,”)
Vis6afp=childfp + str(subj) + '_rest_VIS6A.csv'
Vis6Amp_Agg=np.savetxt(Vis6Amp_Agg, Vis6afp,delimiter=“,”)
Vis10afp=childfp + str(subj) + '_rest_VIS10A.csv'
Vis10Amp_Agg=np.savetxt(Vis10Amp_Agg, Vis10afp,delimiter=“,”)
Aud16afp=childfp + str(subj) + '_rest_AUD16A.csv'
Aud16Amp_Agg=np.savetxt(Aud16Amp_Agg, Aud16afp,delimiter=“,”)

##### save out distribution of DMNpeak-DMNpeak spans
PeakIntervalSpansfp=childfp + str(subj) + '_rest_PeakIntervalSpan.csv'
PeakIntervalSpans=np.savetxt(PeakIntervalSpans, PeakIntervalSpansfp,delimiter=“,”)


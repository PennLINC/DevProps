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
DM_inds12=np.where(ParcelC==12)[1]
DM_inds1=np.where(ParcelC==1)[1]
DM_inds8=np.where(ParcelC==8)[1]
DM_inds=np.hstack((DM_inds12,DM_inds1,DM_inds8))
Mot_inds2=np.where(ParcelC==2)[1]
Mot_inds4=np.where(ParcelC==4)[1]
Mot_inds11=np.where(ParcelC==11)[1]
Mot_inds13=np.where(ParcelC==13)[1]
Vis_inds6=np.where(ParcelC==6)[1]
Vis_inds10=np.where(ParcelC==10)[1]
DAN_inds14=np.where(ParcelC==14)[1]
DAN_inds5=np.where(ParcelC==5)[1]
VAN_inds9=np.where(ParcelC==9)[1]
VAN_inds7=np.where(ParcelC==7)[1]
FPN_inds3=np.where(ParcelC==3)[1]
FPN_inds15=np.where(ParcelC==15)[1]
FPN_inds17=np.where(ParcelC==17)[1]
Aud_inds=np.where(ParcelC==16)[1]

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
	subjDataCor_inSeg=subjDataCort[SegStartInd:(SegStartInd+SegSpan-1),:]
	# get time series assoc. with each parcel
	# MEAN DM
	meanDM=np.mean(subjDataCor_inSeg[:,DM_inds],axis=1)
	DM=meanDM
	# DM
	meanDM1=np.mean(subjDataCor_inSeg[:,DM_inds1],axis=1)
	DM1=meanDM1
	meanDM8=np.mean(subjDataCor_inSeg[:,DM_inds8],axis=1)
	DM8=meanDM8
	meanDM12=np.mean(subjDataCor_inSeg[:,DM_inds12],axis=1)
	DM12=meanDM12
	# FP
	meanFPN3=np.mean(subjDataCor_inSeg[:,FPN_inds3],axis=1)
	FP3=meanFPN3
	meanFPN15=np.mean(subjDataCor_inSeg[:,FPN_inds15],axis=1)
	FP15=meanFPN15
	meanFPN17=np.mean(subjDataCor_inSeg[:,FPN_inds17],axis=1)
	FP17=meanFPN17
	# VAN
	meanVAN7=np.mean(subjDataCor_inSeg[:,VAN_inds7],axis=1)
	VAN7=meanVAN7
	meanVAN9=np.mean(subjDataCor_inSeg[:,VAN_inds9],axis=1)
	VAN9=meanVAN9
	# DAN
	meanDAN5=np.mean(subjDataCor_inSeg[:,DAN_inds5],axis=1)
	DAN5=meanDAN5
	meanDAN14=np.mean(subjDataCor_inSeg[:,DAN_inds14],axis=1)
	DAN14=meanDAN14
	# MOT
	meanMot2=np.mean(subjDataCor_inSeg[:,Mot_inds2],axis=1)
	Mot2=meanMot2
	meanMot4=np.mean(subjDataCor_inSeg[:,Mot_inds4],axis=1)
	Mot4=meanMot4
	meanMot11=np.mean(subjDataCor_inSeg[:,Mot_inds11],axis=1)
	Mot11=meanMot11
	meanMot13=np.mean(subjDataCor_inSeg[:,Mot_inds13],axis=1)
	Mot13=meanMot13
	# VIS
	meanVis6=np.mean(subjDataCor_inSeg[:,Vis_inds6],axis=1)
	Vis6=meanVis6
	meanVis10=np.mean(subjDataCor_inSeg[:,Vis_inds10],axis=1)
	Vis10=meanVis10
	# AUD
	meanAud=np.mean(subjDataCor_inSeg[:,Aud_inds],axis=1)
	AUD16=meanAud
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
	# init amp
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
		DM1Peak=signal.find_peaks(DM1InInterval,distance=IntervalSpan)[0]-1
		DM8Peak=signal.find_peaks(DM8InInterval,distance=IntervalSpan)[0]-1
		DM12Peak=signal.find_peaks(DM12InInterval,distance=IntervalSpan)[0]-1
		FP3Peak=signal.find_peaks(FP3InInterval,distance=IntervalSpan)[0]-1
		FP15Peak=signal.find_peaks(FP15InInterval,distance=IntervalSpan)[0]-1
		FP17Peak=signal.find_peaks(FP17InInterval,distance=IntervalSpan)[0]-1
		VAN7Peak=signal.find_peaks(VAN7InInterval,distance=IntervalSpan)[0]-1
		VAN9Peak=signal.find_peaks(VAN9InInterval,distance=IntervalSpan)[0]-1
		DAN5Peak=signal.find_peaks(DAN5InInterval,distance=IntervalSpan)[0]-1
		DAN14Peak=signal.find_peaks(DAN14InInterval,distance=IntervalSpan)[0]-1
		mot2Peak=signal.find_peaks(Mot2InInterval,distance=IntervalSpan)[0]-1
		mot4Peak=signal.find_peaks(Mot4InInterval,distance=IntervalSpan)[0]-1
		mot11Peak=signal.find_peaks(Mot11InInterval,distance=IntervalSpan)[0]-1
		mot13Peak=signal.find_peaks(Mot13InInterval,distance=IntervalSpan)[0]-1
		vis6Peak=signal.find_peaks(Vis6InInterval,distance=IntervalSpan)[0]-1
		vis10Peak=signal.find_peaks(Vis10InInterval,distance=IntervalSpan)[0]-1
		aud16Peak=signal.find_peaks(Aud16InInterval,distance=IntervalSpan)[0]-1
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
	DM1phases_Agg=np.concatenate((DM1phases_Agg, DM1phases))
	DM8phases_Agg=np.concatenate((DM8phases_Agg, DM8phases))
	DM12phases_Agg=np.concatenate((DM12phases_Agg, DM12phases))
	FP3phases_Agg=np.concatenate((FP3phases_Agg, FP3phases))
	FP15phases_Agg=np.concatenate((FP15phases_Agg, FP15phases))
	FP17phases_Agg=np.concatenate((FP17phases_Agg, FP17phases))
	VAN7phases_Agg=np.concatenate((VAN7phases_Agg, VAN7phases))
	VAN9phases_Agg=np.concatenate((VAN9phases_Agg, VAN9phases))
	DAN5phases_Agg=np.concatenate((DAN5phases_Agg, DAN5phases))
	DAN14phases_Agg=np.concatenate((DAN14phases_Agg, DAN14phases))
	Mot2phases_Agg=np.concatenate((Mot2phases_Agg, Mot2phases))
	Mot4phases_Agg=np.concatenate((Mot4phases_Agg, Mot4phases))
	Mot11phases_Agg=np.concatenate((Mot11phases_Agg, Mot11phases))
	Mot13phases_Agg=np.concatenate((Mot13phases_Agg, Mot13phases))
	Vis6phases_Agg=np.concatenate((Vis6phases_Agg, Vis6phases))
	Vis10phases_Agg=np.concatenate((Vis10phases_Agg, Vis10phases))
	Aud16phases_Agg=np.concatenate((Aud16phases_Agg, Aud16phases))
	# amplitude
	DM1Amp_Agg=np.concatenate((DM1Amp_Agg, DM1Amp))
	DM8Amp_Agg=np.concatenate((DM8Amp_Agg, DM8Amp))
	DM12Amp_Agg=np.concatenate((DM12Amp_Agg, DM12Amp))
	FP3Amp_Agg=np.concatenate((FP3Amp_Agg, FP3Amp))
	FP15Amp_Agg=np.concatenate((FP15Amp_Agg, FP15Amp))
	FP17Amp_Agg=np.concatenate((FP17Amp_Agg, FP17Amp))
	VAN7Amp_Agg=np.concatenate((VAN7Amp_Agg, VAN7Amp))
	VAN9Amp_Agg=np.concatenate((VAN9Amp_Agg, VAN9Amp))
	DAN5Amp_Agg=np.concatenate((DAN5Amp_Agg, DAN5Amp))
	DAN14Amp_Agg=np.concatenate((DAN14Amp_Agg, DAN14Amp))
	Mot2Amp_Agg=np.concatenate((Mot2Amp_Agg, Mot2Amp))
	Mot4Amp_Agg=np.concatenate((Mot4Amp_Agg, Mot4Amp))
	Mot11Amp_Agg=np.concatenate((Mot11Amp_Agg, Mot11Amp))
	Mot13Amp_Agg=np.concatenate((Mot13Amp_Agg, Mot13Amp))
	Vis6Amp_Agg=np.concatenate((Vis6Amp_Agg, Vis6Amp))
	Vis10Amp_Agg=np.concatenate((Vis10Amp_Agg, Vis10Amp))
	Aud16Amp_Agg=np.concatenate((Aud16Amp_Agg, Aud16Amp))
	# same for DM peak span
	PeakIntervalSpans=np.concatenate((PeakIntervalSpans,PeaksandDurs[:,1]))
######## end for each segment
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
##### save out phase and magnitude info
DM1pfp=childfp + str(subj) + '_rest_DM1P.csv'
np.savetxt(DM1pfp, DM1phases_Agg, delimiter=",")
DM8pfp=childfp + str(subj) + '_rest_DM8P.csv'
np.savetxt(DM8pfp, DM8phases_Agg, delimiter=",")
DM12pfp=childfp + str(subj) + '_rest_DM12P.csv'
np.savetxt(DM12pfp, DM12phases_Agg, delimiter=",")
FP3pfp=childfp + str(subj) + '_rest_FP3P.csv'
np.savetxt(FP3pfp, FP3phases_Agg, delimiter=",")
FP15pfp=childfp + str(subj) + '_rest_FP15P.csv'
np.savetxt(FP15pfp, FP15phases_Agg, delimiter=",")
FP17pfp=childfp + str(subj) + '_rest_FP17P.csv'
np.savetxt(FP17pfp, FP17phases_Agg, delimiter=",")
VAN7pfp=childfp + str(subj) + '_rest_VAN7P.csv'
np.savetxt(VAN7pfp, VAN7phases_Agg, delimiter=",")
VAN9pfp=childfp + str(subj) + '_rest_VAN9P.csv'
np.savetxt(VAN9pfp, VAN9phases_Agg, delimiter=",")
DAN5pfp=childfp + str(subj) + '_rest_DAN5P.csv'
np.savetxt(DAN5pfp, DAN5phases_Agg, delimiter=",")
DAN14pfp=childfp + str(subj) + '_rest_DAN14P.csv'
np.savetxt(DAN14pfp, DAN14phases_Agg, delimiter=",")
Mot2pfp=childfp + str(subj) + '_rest_MOT2P.csv'
np.savetxt(Mot2pfp, Mot2phases_Agg, delimiter=",")
Mot4pfp=childfp + str(subj) + '_rest_MOT4P.csv'
np.savetxt(Mot4pfp, Mot4phases_Agg, delimiter=",")
Mot11pfp=childfp + str(subj) + '_rest_MOT11P.csv'
np.savetxt(Mot11pfp, Mot11phases_Agg, delimiter=",")
Mot13pfp=childfp + str(subj) + '_rest_MOT13P.csv'
np.savetxt(Mot13pfp, Mot13phases_Agg, delimiter=",")
Vis6pfp=childfp + str(subj) + '_rest_VIS6P.csv'
np.savetxt(Vis6pfp, Vis6phases_Agg, delimiter=",")
Vis10pfp=childfp + str(subj) + '_rest_VIS10P.csv'
np.savetxt(Vis10pfp, Vis10phases_Agg, delimiter=",")
Aud16pfp=childfp + str(subj) + '_rest_AUD16P.csv'
np.savetxt(Aud16pfp, Aud16phases_Agg, delimiter=",")

# amplitudes
DM1afp=childfp + str(subj) + '_rest_DM1A.csv'
np.savetxt(DM1afp, DM1Amp_Agg, delimiter=",")
DM8afp=childfp + str(subj) + '_rest_DM8A.csv'
np.savetxt(DM8afp, DM8Amp_Agg, delimiter=",")
DM12afp=childfp + str(subj) + '_rest_DM12A.csv'
np.savetxt(DM12afp, DM12Amp_Agg, delimiter=",")
FP3afp=childfp + str(subj) + '_rest_FP3A.csv'
np.savetxt(FP3afp, FP3Amp_Agg, delimiter=",")
FP15afp=childfp + str(subj) + '_rest_FP15A.csv'
np.savetxt(FP15afp, FP15Amp_Agg, delimiter=",")
FP17afp=childfp + str(subj) + '_rest_FP17A.csv'
np.savetxt(FP17afp, FP17Amp_Agg, delimiter=",")
VAN7afp=childfp + str(subj) + '_rest_VAN7A.csv'
np.savetxt(VAN7afp, VAN7Amp_Agg, delimiter=",")
VAN9afp=childfp + str(subj) + '_rest_VAN9A.csv'
np.savetxt(VAN9afp, VAN9Amp_Agg, delimiter=",")
DAN5afp=childfp + str(subj) + '_rest_DAN5A.csv'
np.savetxt(DAN5afp, DAN5Amp_Agg, delimiter=",")
DAN14afp=childfp + str(subj) + '_rest_DAN14A.csv'
np.savetxt(DAN14afp, DAN14Amp_Agg, delimiter=",")
Mot2afp=childfp + str(subj) + '_rest_MOT2A.csv'
np.savetxt(Mot2afp, Mot2Amp_Agg, delimiter=",")
Mot4afp=childfp + str(subj) + '_rest_MOT4A.csv'
np.savetxt(Mot4afp, Mot4Amp_Agg, delimiter=",")
Mot11afp=childfp + str(subj) + '_rest_MOT11A.csv'
np.savetxt(Mot11afp, Mot11Amp_Agg, delimiter=",")
Mot13afp=childfp + str(subj) + '_rest_MOT13A.csv'
np.savetxt(Mot13afp, Mot13Amp_Agg, delimiter=",")
Vis6afp=childfp + str(subj) + '_rest_VIS6A.csv'
np.savetxt(Vis6afp, Vis6Amp_Agg, delimiter=",")
Vis10afp=childfp + str(subj) + '_rest_VIS10A.csv'
np.savetxt(Vis10afp, Vis10Amp_Agg, delimiter=",")
Aud16afp=childfp + str(subj) + '_rest_AUD16A.csv'
np.savetxt(Aud16afp, Aud16Amp_Agg, delimiter=",")

##### save out distribution of DMNpeak-DMNpeak spans
PeakIntervalSpansfp=childfp + str(subj) + '_rest_PeakIntervalSpan.csv'
np.savetxt(PeakIntervalSpansfp, PeakIntervalSpans, delimiter=",")


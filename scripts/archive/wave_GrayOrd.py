import nibabel as nb
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import hcp_utils as hcp

# Subject is set to the 1st passed argument
subj = sys.argv[1]
# wave number is set to the second
waveNum = sys.argv[2]
# parent filepath
parentfp='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' + str(subj) + '/ses-baselineYear1Arm1/func/'
# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'
# load in PG
subjPGfn=childfp+str(subj)+'_PG_LR_32k.dscalar.nii'
PG=nb.load(subjPGfn)
# import waveTRs table (just resting state for now)
waveTabfn=childfp+'/'+str(subj)+'_rest_waveTRs.csv'
waveTab=np.genfromtxt(waveTabfn,delimiter=',')
# import motion/outlier masked AND valid segments only dtseries
filepath=parentfp + str(subj) + '_p2mm_masked_filtered_rest.dtseries.nii'
subjData=nb.load(filepath)
print(int(waveNum))
# match wave requested to row in waveTRs table
wave=waveTab[int(waveNum),:]
# extract specified time interval from time series
start=int(wave[0]+wave[1]-1)
end=int(wave[0]+wave[2]-1)
# and start 3 TRs before and extend 3 TRs after
start=start-3
end=end+3
# sep left and right hemi
PGL=PG.dataobj[0,hcp.struct.cortex_left]
SDL=subjData.dataobj[start:end,hcp.struct.cortex_left]
# and right 1st pg
PGR=PG.dataobj[0,hcp.struct.cortex_right]
SDR=subjData.dataobj[start:end,hcp.struct.cortex_right]
# order time series in accord. w/ pg
# flip so top of hierarchy is top of plot
PGindicesL=np.flip(np.argsort(PGL))
PGindicesR=np.flip(np.argsort(PGR))
# order subj data hierarchically
SDL=SDL[:,PGindicesL]
SDR=SDR[:,PGindicesR]
# convert to percentiles
Lbins=np.zeros(((end-start),70))
Rbins=np.zeros(((end-start),70))
for b in range(70):
	# 1.42857 is modifier to allow "70" to fill 70 percentile bins all the way up to 100
	gradPrctile=np.percentile(PG.dataobj[0,hcp.struct.cortex],(b*1.42857))
	gradPrctile_upper=np.percentile(PG.dataobj[0,hcp.struct.cortex],(b*1.42857)+1.42857)
	# index of vertices belonging to this percentile
	boolean_of_interestL=np.logical_and(PGL > gradPrctile, PGL < gradPrctile_upper)
	boolean_of_interestR=np.logical_and(PGR > gradPrctile, PGR < gradPrctile_upper)
	PGindicesL=np.nonzero(boolean_of_interestL)
	PGindicesR=np.nonzero(boolean_of_interestR)
	PGindices_arrayL=np.array(PGindicesL)
	PGindices_arrayR=np.array(PGindicesR)
	# initialize array of all these vertices to average over
	initPGbinL=SDL[:,PGindices_arrayL]
	initPGbinR=SDR[:,PGindices_arrayR]
	# average signal over this bin 
	meanSigL=np.mean(initPGbinL,axis=2)
	meanSigL=meanSigL[:,0]
	meanSigR=np.mean(initPGbinR,axis=2)
	meanSigR=meanSigR[:,0]
	# plop into ProcTS_bins
	Lbins[:,b]=meanSigL
	Rbins[:,b]=meanSigR
# normalize
Lbins=Lbins-np.mean(Lbins,0)
Lstds=np.mean(np.std(Lbins))
Lbins=Lbins/Lstds
Rbins=Rbins-np.mean(Rbins,0)
Rstds=np.mean(np.std(Rbins))
Rbins=Rbins/Rstds
# transpose
RbinsT=np.transpose(Lbins)
LbinsT=np.transpose(Rbins)
# fig name
fn=childfp+str(subj)+'_rest_wave_gOrd_'+str(waveNum)+'.png'
# make grayOrd plot at the specified time interval
plt.rcParams['figure.dpi'] = 500
plt.subplot(1,2,1)
fig1=plt.matshow(LbinsT,cmap='gray',vmin=-2, vmax=2,extent=[0,(end-start),0,1], aspect=120,fignum=False)
plt.subplot(1,2,2)
fig2=plt.matshow(RbinsT,cmap='gray',vmin=-2, vmax=2,extent=[0,(end-start),0,1], aspect=120,fignum=False)
plt.savefig(fn,bbox_inches='tight')

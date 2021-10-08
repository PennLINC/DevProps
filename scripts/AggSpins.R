# load in subjects and arguments
Args=commandArgs(trailingOnly=TRUE)
subj=Args[1]
SpinNum=Args[2]
paste(subj)
paste(SpinNum)
# set child filepath
cfp=paste0('/cbica/projects/abcdfnets/results/wave_output/',subj,'/')
# set tasks
tasks=c('rest','SST','nback','MID')
for (t in 1:4){
	task=tasks[t]
	# load in pre-intialized matrices
	MagPGbinFn=paste0(cfp,task,'MagPGBinSpin.csv')
	PhasePGbinFn=paste0(cfp,task,'PhasePGBinSpin.csv')
	DurationFn=paste0(cfp,task,'AvgDurSpin.csv')
	MagPGbin=read.csv(MagPGbinFn,header=F)
	PhasePGbin=read.csv(PhasePGbinFn,header=F)
	Duration=read.csv(DurationFn,header=F)
	# load in metrics calculated for these subjs
	WMagFn=paste0(cfp,subj,'_',tasks[t],'_MagMat_Thr_spun.csv')
	WDelayFn=paste0(cfp,subj,'_',tasks[t],'_delayMat_Thr_spun.csv')
	WDurFn=paste0(cfp,subj,'_',tasks[t],'_waveTRs_spun.csv')
	WMag=read.csv(WMagFn,header=F)
	DelayMat=read.csv(WDelayFn,header=F)
	waveTRs=read.csv(WDurFn,header=F)
	# convert Magnitude to average values per pg bin
	WMag[is.na(WMag)]=0
	BinWiseMag=rowMeans(WMag)
	MagPGbin[SpinNum,]=t(BinWiseMag)
	# save updated magnitude back out
	write.table(MagPGbin,MagPGbinFn,row.names=F,col.names=F)
	### Done with magnitude
	# make a column representing when the top-o-pg peak is b/w its troughs
	waveTRs$TPGpeak=waveTRs$V2+waveTRs$V4
	# make a column representing the lower bound of each cycle
	waveTRs$LB=waveTRs$V2-waveTRs$TPGpeak
	# same for upper bound
	waveTRs$UB=waveTRs$V3-waveTRs$TPGpeak
	# determine wave timepoints from duration file
	startEnddif=waveTRs$V3-waveTRs$V2
	# get median speed of wave
	MedSped=median(startEnddif)
	# write out median speed
	Duration[SpinNum]=MedSped
	write.table(Duration,DurationFn,row.names=F,col.names=F)
	# CAN ADD TIME SPENT IN WAVE FOR NULL COMPARISONS IF NEEDED
	# initialize phase-PGBin matrix - number of waves + number of non-top bins
	PhPGB=matrix(999,length(startEnddif),24)
	# average/median for phases, need to delineate on wave by wave basis
	for (W in 1:length(startEnddif)){
		# just removing self-referential 25th top-o-pg bin
		DelayForThisWave=DelayMat[1:24,W]
		# negatives need to be relative to lower bound
		NegInds=which(DelayForThisWave<0)
		# get how far (percentagewise) this peak precedes in top-o-pg cycle, retain negativity for parsing later
		PhPGB[W,NegInds]=(DelayForThisWave[NegInds])/(waveTRs$LB[W])*-1
		# positives need to be relative to upper bound
		PosInds=which(DelayForThisWave>0)
		PhPGB[W,PosInds]=(DelayForThisWave[PosInds])/(waveTRs$UB[W])
		# note those in sync
		synced=which(DelayForThisWave==0)
		PhPGB[W,synced]=0
		# note NaNs
		NaInds=is.na(DelayForThisWave)
		PhPGB[W,NaInds]=NaN
	}
	# extract mean phase offset for each PGBin
	phOff=colMeans(abs(PhPGB),na.rm=T)
	# put back in
	PhasePGbin[SpinNum,]=t(phOff)
	# save back out
	write.table(PhasePGbin,PhasePGbinFn,row.names=F,col.names=F)
	# and just average duration
}

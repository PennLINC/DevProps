subj=commandArgs(trailingOnly=TRUE)
# for each task
tasks=c('rest','SST','nback','MID')
# common output Dir
childfp=paste0('/cbica/projects/abcdfnets/results/wave_output/',subj,'/')
# for each task
for (t in 1:4){
	# Mag from magMat
	MagMatFN=paste0(childfp,subj,'_',tasks[t],'_MagMat_Thr.csv')
	mtry <- try(read.csv(MagMatFN, header = F),silent=T)
	if (class(mtry) == "try-error") {
		message(paste(MagMatFN,' is empty'))
	} else { MagMat=read.csv(MagMatFN,header=F)
		# derive binWise Mag
		MagMat[is.na(MagMat)]=0
		BinWiseMag=rowMeans(MagMat)
		# save binWise Mag
		write.table(BinWiseMag,paste0(childfp,subj,'_',tasks[t],'_BinWisePWMagnitude.csv'),row.names=F,col.names=F)
		# save normative Mag
		write.table(mean(BinWiseMag),paste0(childfp,subj,'_',tasks[t],'_MeanPWMagnitude.csv'),row.names=F,col.names=F)
	}
	# get speed and time spent in waves from waveTRs; phase from waveTRs and relative delay
	waveTRsFN=paste0(childfp,subj,'_',tasks[t],'_waveTRs.csv')
	DelayMatFN=paste0(childfp,subj,'_',tasks[t],'_delayMat_Thr.csv')
	mtry <- try(read.csv(waveTRsFN, header = F),silent=T)
        if (class(mtry) == "try-error") {
                message(paste(waveTRsFN,' is empty'))
        } else { waveTRs=read.csv(waveTRsFN,header=F)
                # derive time spent in waves
		startEnddif=waveTRs$V3-waveTRs$V2
		# Total time spent in waves, total surviving TRs saved elsewhere in pipeline for comparison
                TTSW=sum(startEnddif)
		# save total time
                write.table(TTSW,paste0(childfp,subj,'_',tasks[t],'_TotPWTime.csv'),row.names=F,col.names=F)
		# Speed
		MedSped=median(startEnddif)
		# save median speed
		write.table(MedSped,paste0(childfp,subj,'_',tasks[t],'_MedPWDuration.csv'),row.names=F,col.names=F)
		# load delay matrix, should exist if waveTRs does
		DelayMat=read.csv(DelayMatFN,header=F)
		# make a column representing when the top-o-pg peak is b/w its troughs
		waveTRs$TPGpeak=waveTRs$V2+waveTRs$V4
		# make a column representing the lower bound of each cycle
		waveTRs$LB=waveTRs$V2-waveTRs$TPGpeak
		# same for upper bound
		waveTRs$UB=waveTRs$V3-waveTRs$TPGpeak
		# initialize phase-PGBin matrix - number of waves + number of non-top bins
		PhPGB=matrix(999,length(startEnddif),24)
		# for each wave, get % across duration that peak occurs w/in each bin
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
			# there should be no 999s left
        	}
		write.table(PhPGB,paste0(childfp,subj,'_',tasks[t],'_BinWisePhase.csv'),row.names=F,col.names=F)
	}
}

### aggregate subject phase by pg bin
# subject list
subjs=read.delim('~/clean100.txt',header=F)
# for each task
tasks=c('rest','SST','nback','MID')
# common output dir
childfp=paste0('/cbica/projects/abcdfnets/results/wave_output/')
# for each task
for (t in 1:4){
	# initialize output df, extra col for subj name
	phaseNbins=data.frame(matrix(0,length(subjs$V1),25))
	# populate with each subj
	for (s in 1:length(subjs$V1)){
		# subjects s
		subj=subjs$V1[s]
		phaseNbins[s,1]=subj	
		# load phase by PG bin
		phBpgFN=paste0(childfp,subj,'/',subj,'_',tasks[t],'_BinWisePhase.csv')
		phBpg=read.delim(phBpgFN,header=F,sep=' ')
		# extract mean phase offset from each PGbin
		phOff=colMeans(abs(phBpg),na.rm=T)
		phaseNbins[s,2:25]=phOff
	}
	# save it
	write.csv(phaseNbins,paste0('/cbica/projects/abcdfnets/results/',tasks[t],'_BinWiseMeanPhase.csv'))
}

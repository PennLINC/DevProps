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
	MagNbins=data.frame(matrix(0,length(subjs$V1),26))
	# populate with each subj
	for (s in 1:length(subjs$V1)){
		# subjects s
		subj=subjs$V1[s]
		MagNbins[s,1]=subj	
		# load Mag by PG bin
		MBpgFN=paste0(childfp,subj,'/',subj,'_',tasks[t],'_BinWisePWMagnitude.csv')
		MBpg=read.delim(MBpgFN,header=F,sep=' ')
		# extract mean Mag from each PGbin
		MMag=MBpg$V1
		MagNbins[s,2:26]=MMag
	}
	# save it
	write.csv(MagNbins,paste0('/cbica/projects/abcdfnets/results/',tasks[t],'_BinWiseMeanMag.csv'))
}

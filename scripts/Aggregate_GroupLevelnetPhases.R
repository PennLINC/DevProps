### Aggregate group-level netPhase stats
# output directory to pull from
Odir='/cbica/projects/abcdfnets/results/wave_output'
# set names of all metrics being aggregated
NetMets<-c('DM1P','MOT2P','FP3P','MOT4P','DAN5P','VIS6P','VAN7P','DM8P','VAN9P','VIS10P','MOT11P','DM12P','MOT13P','DAN14P','FP15P','AUD16P','FP17P','DM1A','MOT2A','FP3A','MOT4A','DAN5A','VIS6A','VAN7A','DM8A','VAN9A','VIS10A','MOT11A','DM12A','MOT13A','DAN14A','FP15A','AUD16A','FP17A','PeakIntervalSpan')
# get subjects ran in this iteration
subjs=read.delim('~/clean100.txt',header=F)
# init list
theList <- list()
# for each metric
for (m in 1:length(NetMets)){
	# get name of this metrics
	Metric=NetMets[m]
	MetricList <- list()
	# load every subject in
	for (s in 1:length(subjs$V1)){
		filename=paste0(Odir,'/',subjs$V1[s],'/',subjs$V1[s],'_rest_',Metric,'.csv')
		file=read.csv(filename,header=F)
		# port into list
		MetricList <- append(MetricList,file$V1)
	}
	theList[[m]]=MetricList
}
df<-do.call("cbind",theList)
colnames(df)<-NetMets
# save out group-level arrays
write.csv(df,'/cbica/projects/abcdfnets/results/restGroup-level_netphase.csv',row.names=F)

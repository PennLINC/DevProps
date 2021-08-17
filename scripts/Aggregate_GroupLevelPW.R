### Aggregate group-level principal wave stats

# define a function to load in each subj file and stack wave instances on prev. subjs
#LoadNstack <- function(x){
	# for each listed file
#	file=paste0(x,'/',list.files(x,pattern=fileExtn))
#	loadedFile=read.csv(file,header=F)
        # transpose loaded file so each row is a wave
#        loadedFile=t(loadedFile)
#        colnames(loadedFile)=c(colnames(InitDf),'TRs')
        # make next W rows this subject's data
#        InitDf=rbind(InitDf,loadedFile[,1:4])
#}

# for second PG
#LoadNstack2 <- function(x){
        # for each listed file
#        file=paste0(x,'/',list.files(x,pattern=fileExtn2))
#        loadedFile=read.csv(file,header=F)
        # transpose loaded file so each row is a wave
#	loadedFile=t(loadedFile)
#	colnames(loadedFile)=c(colnames(InitDf2),'TRs')
	# make next W rows this subject's data
#        InitDf2=rbind(InitDf2,loadedFile[,1:4])
#}



# For each task
tasks = c('rest','sst','nback','mid')

for (i in 1:4){
	# initialize dataframeto tack a new row onto for each wave instance
	InitDf<-data.frame(0,0,0,0)	
	InitDf2<-data.frame(0,0,0,0)  
	# column names to match output structure
	colnames(InitDf)<-c('PGCor','TempSpan','EarliestPG_Bin','RelMagSlope')
	colnames(InitDf2)<-c('PG2Cor','TempSpan2','EarliestPG2_Bin','RelMagSlope2')
	# file extension of this task
	fileExtn=paste0('*',tasks[i],'_waveProps.csv')
	fileExtn2=paste0('*',tasks[i],'_waveProps_PG2.csv')
	# get filenames of this task for pg1
	dirs=list.dirs('/cbica/projects/abcdfnets/results/wave_output')
	# remove first dir, is parent dir
	dirs=dirs[-1]
	# take this tasks' wave recordings out of directories, stack	
	for (s in 1:length(dirs)){
		if (length(list.files(dirs[s],pattern=fileExtn)) != 0) {
	        	file=paste0(dirs[s],'/',list.files(dirs[s],pattern=fileExtn))
			loadedFile=read.csv(file,header=F)
        		# transpose loaded file so each row is a wave
        		loadedFile=t(loadedFile)
        		colnames(loadedFile)=c(colnames(InitDf),'TRs')
        		# make next W rows this subject's data
        		InitDf=rbind(InitDf,loadedFile[,1:4])
			}
	}
	#pg2
        for (s in 1:length(dirs)){
                if (length(list.files(dirs[s],pattern=fileExtn2)) != 0) {	
			file=paste0(dirs[s],'/',list.files(dirs[s],pattern=fileExtn2))
			loadedFile=read.csv(file,header=F)
               		# transpose loaded file so each row is a wave
                	loadedFile=t(loadedFile)
                	colnames(loadedFile)=c(colnames(InitDf2),'TRs')
                	# make next W rows this subject's data
                	InitDf2=rbind(InitDf2,loadedFile[,1:4])
        	}
	}
	# save out group-level arrays
	saveRDS(InitDf,paste0('/cbica/projects/abcdfnets/results/',tasks[i],'Group-level_waves'))
	saveRDS(InitDf2,paste0('/cbica/projects/abcdfnets/results/',tasks[i],'Group-level_waves2'))  	
}

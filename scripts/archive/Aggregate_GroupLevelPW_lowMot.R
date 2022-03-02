### Aggregate group-level principal wave stats

# For each task
tasks = c('rest','SST','nback','MID')

# read in low motion subjects
LowMotSubjs=readRDS('/cbica/projects/abcdfnets/dropbox/vLowMotSubjs.rds')

# individualized PG
for (i in 1:4){
	# initialize dataframeto tack a new row onto for each wave instance
	InitDf<-data.frame(0,0,0,0)	
	InitDf2<-data.frame(0,0,0,0)  
	# column names to match output structure
	colnames(InitDf)<-c('PGCor','TempSpan','EarliestPG_Bin','RelMagSlope')
	colnames(InitDf2)<-c('PGCor','TempSpan','EarliestPG_Bin','RelMagSlope')
	# file extension of this task
	fileExtn=paste0('*',tasks[i],'_waveProps.csv')
	fileExtn2=paste0('*',tasks[i],'_waveProps_gPG.csv')
	# get filenames of this task for pg1
	dirs=list.dirs('/cbica/projects/abcdfnets/results/wave_output')
	# remove first dir, is parent dir
	dirs=dirs[-1]
	# take this tasks' wave recordings out of directories, stack	
	for (s in 1:length(dirs)){
		if (length(list.files(dirs[s],pattern=fileExtn)) != 0) {
	        	# adding contingency to only read if subject in lowMot select
			strings=strsplit(dirs[s],'/')
			subj=strings[[1]][7]
			if (subj %in% LowMotSubjs) {
				message(paste0(subj,' low motion subj'))
				file=paste0(dirs[s],'/',list.files(dirs[s],pattern=fileExtn))
				# try to read, record 0-line instances for subjects with no PWs
				mtry <- try(read.csv(file, header = F),silent=T)
				if (class(mtry) != "try-error") { 
					loadedFile=read.csv(file, header = F)
        				# transpose loaded file so each row is a wave
        				loadedFile=t(loadedFile)
        				colnames(loadedFile)=c(colnames(InitDf),'TRs')
        				# make next W rows this subject's data
        				InitDf=rbind(InitDf,loadedFile[,1:4])
				} else {
					message(paste(file, ": 0-line file"))	
				}
			} else {
				message(paste0(subj,' dont meet no motion threshold'))
			}
		}
	}

# group PG
        for (s in 1:length(dirs)){
                if (length(list.files(dirs[s],pattern=fileExtn2)) != 0) {	
			# adding contingency to only read if subject in lowMot select
                        strings=strsplit(dirs[s],'/')
                        subj=strings[[1]][7]
                        if (subj %in% LowMotSubjs) {
				message(paste0(subj,' low motion subj'))
				file=paste0(dirs[s],'/',list.files(dirs[s],pattern=fileExtn2))
				# try as used above
				mtry2 <- try(read.csv(file, header = F),silent=T)
				if (class(mtry2) != "try-error") {
					loadedFile=read.csv(file,header=F)
              				# transpose loaded file so each row is a wave
                			loadedFile=t(loadedFile)
                			colnames(loadedFile)=c(colnames(InitDf2),'TRs')
                			# make next W rows this subject's data
                			InitDf2=rbind(InitDf2,loadedFile[,1:4])
        			} else {
					message(paste(file, ": 0-line file"))
				}
			} else {
				message(paste0(subj,' dont meet no motion threshold'))
			}
		}
	}
	# save out group-level arrays
	saveRDS(InitDf,paste0('/cbica/projects/abcdfnets/results/',tasks[i],'Group-level_waves_LM'))
	saveRDS(InitDf2,paste0('/cbica/projects/abcdfnets/results/',tasks[i],'Group-level_waves2_LM'))  	

}

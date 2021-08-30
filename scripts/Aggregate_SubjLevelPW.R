### Aggregate subject-level principal wave stats

# record num waves, num TD, num BU, num TRs pulled from, avg. wave speed, sd wave speed

# get directory name of each subject
dirs=list.dirs('/cbica/projects/abcdfnets/results/wave_output')
# remove first dir, is parent dir
dirs=dirs[-1]

# initialize dataframe - subject ID, TRs, num waves, num TD, num BU, avg. speed, sd wave speed x 2 (rest and nback)
df<-data.frame(matrix(nrow=length(dirs),ncol=25))
colnames(df)<-c('ID','TRs_r','numW_r','numTD_r','numBU_r','avgS_r','sdS_r','TRs_n','numW_n','numTD_n','numBU_n','avgS_n','sdS_n','TRs_m','numW_m','numTD_m','numBU_m','avgS_m','sdS_m','TRs_s','numW_s','numTD_s','numBU_s','avgS_s','sdS_s')

# column names of wave counts by task for 0-wave instance recording
cnamesWcs=c('numW_r','numW_n','numW_m','numW_s')

# scan types
tasks = c('rest','nback','mid','SST')

# starting colnames
starting_colnames<-c('PGCor','TempSpan','EarliestPG_Bin','RelMagSlope')

# working threshold for bottom up and top down
TDthresh=-.4
BUthresh=.4

df$ID=dirs

# for each subject
for (s in 1:length(dirs)){
	# For each task
	for (i in 1:4){
		fileExtn=paste0('*',tasks[i],'_waveProps.csv')
		# if file exists
		if (length(list.files(dirs[s],pattern=fileExtn)) != 0) {
			file=paste0(dirs[s],'/',list.files(dirs[s],pattern=fileExtn))
		 	# try to read, record 0-line instances for subjects with no PWs
		        mtry <- try(read.csv(file, header = F),silent=T)
			if (class(mtry) == "try-error") {
				# pertinent colname
				pcn=cnamesWcs[i]
				df[s,pcn]=0
				message(paste(file,' is empty'))	
			} else { loadedFile=read.csv(file,header=F)
			# remove columns NA on PG*delay cor or TR span
			loadedFile=data.frame(loadedFile[,colSums(is.na(loadedFile[1:2,]))==0])
			# if i=1 (resting)
			if (i == 1){
				# record used TRs (after motion mask)
				df$TRs_r[s]=loadedFile[5,1]
				# get number of observed waves
				df$numW_r[s]=length(loadedFile)
				# get number of top downs
				df$numTD_r[s]=length(loadedFile[,loadedFile[1,]<TDthresh])
				# get number of bottom ups
				df$numBU_r[s]=length(loadedFile[,loadedFile[1,]>BUthresh])
				# get average speed in TRs
				df$avgS_r[s]=mean(as.numeric(loadedFile[2,]))
				# get speed SD
				df$sdS_r[s]=sd(as.numeric(loadedFile[2,]))
			# if i=2 (nback)
			} else if (i == 2) {
				# record used TRs (after motion mask)
        df$TRs_n[s]=loadedFile[5,1]
        # get number of observed waves
        df$numW_n[s]=length(loadedFile)
        # get number of top downs
        df$numTD_n[s]=length(loadedFile[,loadedFile[1,]<TDthresh])
        # get number of bottom ups
        df$numBU_n[s]=length(loadedFile[,loadedFile[1,]>BUthresh])
        # get average speed in TRs
        df$avgS_n[s]=mean(as.numeric(loadedFile[2,]))
        # get speed SD
        df$sdS_n[s]=sd(as.numeric(loadedFile[2,]))
      # if i=3 (mid)
      } else if (i == 3) {
        # record used TRs (after motion mask)
        df$TRs_m[s]=loadedFile[5,1]
        # get number of observed waves
        df$numW_m[s]=length(loadedFile)
        # get number of top downs
        df$numTD_m[s]=length(loadedFile[,loadedFile[1,]<TDthresh])
        # get number of bottom ups
        df$numBU_m[s]=length(loadedFile[,loadedFile[1,]>BUthresh])
        # get average speed in TRs
        df$avgS_m[s]=mean(as.numeric(loadedFile[2,]))
        # get speed SD
        df$sdS_m[s]=sd(as.numeric(loadedFile[2,]))  
       # if i=4 (sst)
       } else if (i == 4) {
        # record used TRs (after motion mask)
        df$TRs_s[s]=loadedFile[5,1]
        # get number of observed waves
        df$numW_s[s]=length(loadedFile)
        # get number of top downs
        df$numTD_s[s]=length(loadedFile[,loadedFile[1,]<TDthresh])
        # get number of bottom ups
        df$numBU_s[s]=length(loadedFile[,loadedFile[1,]>BUthresh])
        # get average speed in TRs
        df$avgS_s[s]=mean(as.numeric(loadedFile[2,]))
        # get speed SD
        df$sdS_s[s]=sd(as.numeric(loadedFile[2,]))
        }  			
	}
	}	
}
}

# save out group-level arrays
saveRDS(df,paste0('/cbica/projects/abcdfnets/results/Subject-level_waves'))

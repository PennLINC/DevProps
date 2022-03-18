# extract carit task performance from individ. subject files in /cbica/rbc
# note you will need access to that directory to run this, ran from pinesa user rather than project user for that reason


# get list of subjects to get task perf. for
cmatch_subjs=read.table('/cbica/projects/pinesParcels/PWs/rs_cmatch_subs.csv')

# number of subjects
numsubjs=dim(cmatch_subjs)[1]

# initialize df
taskperf=data.frame(cmatch_subjs$x)
colnames(taskperf)<-'SubjID'

# convert to subj ids of directory structure
taskperf$SubjID<-gsub('sub-','HCD',taskperf$SubjID)
taskperf$perf<-rep(2,numsubjs)
taskperf$rt<-rep(0,numsubjs)

# for each subj
for (s in 1:numsubjs){
  subject=taskperf$SubjID[s]
  # get stats filepath
  statsfp=paste0('/cbica/projects/RBC/testing/hcpd/HCPD/',subject,'/',subject,'_V1_MR/unprocessed/tfMRI_CARIT_AP/LINKED_DATA/PSYCHOPY/CARIT_',subject,'_V1_A_run2_stats.csv')
  if (file.exists(statsfp)) {
    # pull out _stats.csv
    stats=read.csv(statsfp)
    # pull out Hit_Prop, plop into df
    taskperf$perf[s]<-stats$Hit_Prop
    taskperf$rt[s]<-stats$Hit_RT
  }
}

# find instances of "2", where file does not exist
no2s <- subset(taskperf, perf != 2) 

# find instances of "na", where file DOES exist, but is populated with 'NA
# appears to be redundant, but in for sanity check
noNans <- subset(no2s, perf != 'NA')

# convert accuracy to performance
trueperf<-noNans$perf/noNans$rt

cor.test(trueperf,noNans$perf)

noNans$perf<-trueperf

# z-score
noNans$perf<-scale(noNans$perf)

# convert subjIDs back
noNans$SubjID<-gsub('HCD','sub-',noNans$SubjID)

# save out
write.csv(noNans,'/cbica/projects/pinesParcels/PWs/hcpd_caritTaskPerf.csv')

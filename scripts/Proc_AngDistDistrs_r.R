.libPaths('/cbica/home/pinesa/Rlibs')
subj=commandArgs(trailingOnly=TRUE)
print(subj)

library(diptest)

# read in ang distrs for this subj
#gangdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_L.csv'))
#gangdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_gAngDist_Masked4_R.csv'))

# process full-brain-merged distributions
#gmeltedL<-unlist(data.frame(gangdL))
#gmeltedR<-unlist(data.frame(gangdR))
#gmelted=c(gmeltedL,gmeltedR)

# dip test statistic
#gDipT=dip.test(gmelted)
#print(gDipT)
#saveRDS(gDipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_gDipTest.rds'))

# repeat with tertile-calculated
gangdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_tertAngDist_Masked4_L.csv'))
gangdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_tertAngDist_Masked4_R.csv'))

# process full-brain-merged distributions
gmeltedL<-unlist(data.frame(gangdL))
gmeltedR<-unlist(data.frame(gangdR))
gmelted=c(gmeltedL,gmeltedR)

# dip test statistic
gDipT=dip.test(gmelted)
print(gDipT)
saveRDS(gDipT,paste0('~/results/PWs/Proced/',subj,'/',subj,'_tertDipTest.rds'))


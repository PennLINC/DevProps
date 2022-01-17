# aggregate angular distributions, save mean angular distance per subj

# read subject list
subjList=read.delim('~/PWs/hcpd_subj_list.txt')

for (s in 78:dim(subjList)[1]){
  subj=subjList[s,1]
  print(s)
  # if this subject ran, add it to the shape matrix 
  if(file.exists(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))){
 	 # load in data
 	 angdL=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_L.csv'))
 	 angdR=read.csv(paste0('~/results/PWs/Proced/',subj,'/',subj,'_AngDist_Masked_R.csv'))
 	 # melt
 	 meltedL<-unlist(data.frame(angdL))
 	 meltedR<-unlist(data.frame(angdR))
 	 melted=c(meltedL,meltedR)
	 meanAngDist=mean(melted)
	 print(meanAngDist)
	 saveRDS(meanAngDist,(paste0('~/results/PWs/Proced/',subj,'/',subj,'_MeanAngDist')))
}
}

Group_level
================

``` r
library(ggplot2)
library(mgcv)
```

    ## Loading required package: nlme

    ## This is mgcv 1.8-38. For overview type 'help("mgcv-package")'.

``` r
library(scales)
```

``` r
# load in subj list
subjList=read.delim('Y:/PWs/hcpd_subj_list.txt')

# load in ages
demo=read.csv('Y:/PWs/hcpd_demographics.csv')
# convert to used naming convention
demo$SubjID<-gsub('HCD','sub-',demo$src_subject_id)

# load in FD
FD_TRs=read.csv('Y:/PWs/Subj_FD_RemTRs.csv')
colnames(FD_TRs)[1]<-'SubjID'
colnames(FD_TRs)[2]<-'FD'
colnames(FD_TRs)[3]<-'RemainingTRs'

# merge by subjID
mergeddf<-merge(demo,FD_TRs,by='SubjID')

# initialize D
Ds=rep(0,dim(subjList)[1])
Ps=rep(0,dim(subjList)[1])
remTRs=rep(0,dim(subjList)[1])

# subjvec to run in parallel for even more confidence in merging
Subjvec=rep(0,dim(subjList)[1])


# load in D iteratively
for (s in 1:dim(subjList)[1]){
    subj=subjList[s,1]
    ResFp=paste0('Y:/results/PWs/Proced/',subj,'/')
    # if output exists
    if (file.exists(paste0(ResFp,subj,'_gDipTest.rds'))) {
      # load in dip test
      dipstat=readRDS(paste0(ResFp,subj,'_gDipTest.rds'))
      Ds[s]=dipstat$statistic
      Ps[s]=dipstat$p.value
      # and to record same order of subj loading and prevent funny biz
      Subjvec[s]=subj
      # and remaining TRs
      remTRs[s]=mergeddf$RemainingTRs[mergeddf$SubjID==subj]
    }
}

# merge populated vecs
popDf=data.frame(Subjvec,Ds,remTRs)
colnames(popDf)[1]<-'SubjID'

# remaining TRs thresh
inclusionVec<-popDf$remTRs>600
includedSubjsvec<-popDf$SubjID[inclusionVec]
popDfThresh=popDf[popDf$remTRs>600,]

masterdf=merge(mergeddf,popDfThresh,by='SubjID')

# load in dmn seg
DMNseg=read.csv('Y:/results/DMNseg.csv')
colnames(DMNseg)<-c('SubjID','DMNseg','DMNsegGro')
masterdf<-merge(masterdf,DMNseg,by='SubjID')
```

``` r
####### This chunk is just for saving out subject lists corresponding to age tertiles for the cluster

# make list of old subjs vs younger
low_tertile_age=quantile(masterdf$interview_age,.33)
high_tertile_age=quantile(masterdf$interview_age,.66)

# youngins
youngins=masterdf$SubjID[masterdf$interview_age<low_tertile_age]
write.table(youngins,'Y:/PWs/young_subs.csv')
#geezers
geezers=masterdf$SubjID[masterdf$interview_age>high_tertile_age]
write.table(geezers,'Y:/PWs/old_subs.csv')



##### make subject lists for RS w/ good _carit task runs ( > 300 TRs retained)

###### 388
rs_subs<-masterdf$SubjID
#print
write.table(rs_subs,'Y:/PWs/rs_subs.csv')

#### RS_C matched: rope in a bunch of goblyguk just for equivalence with facewise_t
# load in FD
FD_TRs=read.csv('Y:/PWs/Subj_FD_RemTRs.csv')
colnames(FD_TRs)[1]<-'SubjID'
colnames(FD_TRs)[2]<-'FD'
colnames(FD_TRs)[3]<-'RemainingTRs'
# FD for carit runs
FD_TRs_c=read.csv('Y:/PWs/Subj_FD_RemTRs_c.csv')
colnames(FD_TRs_c)[1]<-'SubjID'
colnames(FD_TRs_c)[2]<-'FDc'
colnames(FD_TRs_c)[3]<-'RemainingTRsc'
# merge by subjID
mergeddf<-merge(demo,FD_TRs,by='SubjID')
mergeddf<-merge(mergeddf,FD_TRs_c,by='SubjID')
# exclude subjects with less than 600 TRs remaining (in rest)
inclusionVec<-mergeddf$RemainingTRs>600
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
mergeddf<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(mergeddf)[1]
print(remainingSubjs)
```

    ## [1] 387

``` r
# exclude subjects with less than 300 TRs remaining (in carit)
inclusionVec<-mergeddf$RemainingTRsc>300
# include NAs in the exclusion
inclusionVec[is.na(inclusionVec)==TRUE]=FALSE
# subset the master df accordingly
df<-mergeddf[inclusionVec,]
# get count of remaining subjects
remainingSubjs=dim(df)[1]
print(remainingSubjs)
```

    ## [1] 281

``` r
# extract subjs
rs_cmatch_subjs<-df$SubjID
# print
write.table(rs_cmatch_subjs,'Y:/PWs/rs_cmatch_subs.csv')

# _C: will need FP alteration in .m, but is from same subjects
```

``` r
# read in discretized cross subj histograms (linear)
o_disthist<-t(read.csv('Y:/results/PWs/old_subs_AngDistHist.csv'))
y_disthist<-t(read.csv('Y:/results/PWs/young_subs_AngDistHist.csv'))

# normalize to total numbers
o_disthist<-o_disthist/sum(o_disthist)
y_disthist<-y_disthist/sum(y_disthist)

difHist=o_disthist-y_disthist

# boring plots
#barplot(t(o_disthist),ylim=c(0,0.07),main = 'Old r.s. OpFl Directionality',xlab='0-180')
#barplot(t(y_disthist),ylim=c(0,0.07),main = 'Young OpFl Directionality',xlab='0-180')
#barplot(t(difHist),main='Difference in OpFl Directionality',xlab='0-180')

# fancy plots
#difference fancy

plotdataf<-data.frame(seq(1,18),o_disthist-y_disthist)
# correct to be degrees rather than histogram bin index
plotdataf$seq.1..18.<-plotdataf$seq.1..18.*10

p<-ggplot(plotdataf,aes(x=seq.1..18.,y=o_disthist...y_disthist))
p+geom_bar(stat='identity',color='black',aes(fill= ..x..),binwidth = 10)+xlab("Distance (Degrees)")+theme_classic(base_size = 23)+scale_x_continuous(limits=(c(0,190)))+scale_fill_gradientn("value",colors=c("blue","cyan","green","yellow","orange","red"))+ggtitle("Angular Distance over all TRs: Dif")+theme(plot.title= element_text(size=30, face="bold"), axis.title = element_text(size=30, face="bold"),axis.text = element_text(face="bold",size=30),legend.title=element_blank(),legend.text=element_text(size=20),legend.position=c(1.07,.41),plot.margin=margin(b=.1,t=.1,l=.1,r=2.3, unit='cm'))+ylab('count')
```

    ## Warning: Ignoring unknown parameters: binwidth

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# load in permuted old young dif
PermDif=read.csv('Y:/results/PWs/Old_v_young_AngDistDif_permuted.csv')
permDf<-data.frame(PermDif)

# make vectors of mean values of permuted at each bin
meanVec=colMeans(permDf[1:1000,])
# make vectors of 95 and 5th percentile of permuted at each bin
topPercVec=unlist(lapply(permDf[1:1000,1:18],quantile,prob=.95,names=F))
botPercVec=unlist(lapply(permDf[1:1000,1:18],quantile,prob=.05,names=F))

plotdataf$meanVec=meanVec
plotdataf$topPercVec=topPercVec
plotdataf$botPercVec=botPercVec

plotdf<-data.frame(meanVec,topPercVec,botPercVec,seq(1:18),t(PermDif[1001,]))

p<-ggplot(plotdataf,aes(x=seq.1..18.,y=o_disthist...y_disthist*100))

p+geom_bar(stat='identity',color='black',aes(fill= ..x..),binwidth = 10)+xlab("Distance (Degrees)")+theme_classic(base_size = 26)+scale_x_continuous(limits=(c(5,185)))+scale_fill_gradientn("value",colors=c("blue","cyan","green","yellow","orange","red"))+ggtitle("\u2207PG Distance: Difference")+ylab('Percentage of Propagations')+theme(plot.title= element_text(size=30, face="bold"), axis.title = element_text(size=30, face="bold",vjust=-1),axis.text = element_text(face="bold",size=30),legend.title=element_blank(),legend.text=element_text(size=20),legend.position='none',plot.margin=margin(b=.1,t=.1,l=.1,r=2.3, unit='cm'))+ylab('count')+ylab('Percentage of Propagations')+geom_ribbon(aes(x=seq.1..18.,ymin=botPercVec*100,ymax=topPercVec*100),fill='gray80',alpha=.5)
```

    ## Warning: Ignoring unknown parameters: binwidth

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# saved out at 850x825
```

``` r
# checkup on how BuProp from hists (and rs_subs.csv) maps onto bimodality
BuProp=t(read.csv('Y:/results/PWs/rs_subsBuProp.csv'))
Subjs=read.delim('Y:/PWs/rs_subs.csv',sep=' ')
BuPropdf<-data.frame(BuProp,Subjs)
colnames(BuPropdf)<-c('BuProp','SubjID')
BuMerge<-merge(BuPropdf,masterdf,by='SubjID')

#plotting
BuDs<-ggplot(BuMerge,aes(x=Ds,y=1-BuProp))+geom_point(size=3,alpha=.5)+theme_classic(base_size=28)+geom_smooth(method='lm',color='black')
BuDs+xlab('Dip Statistics')+ylab('Top-Down')+scale_y_continuous(labels = percent_format(accuracy = 1))+scale_x_continuous(breaks=c(.002,.006,.01))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# convert to prop TD
BuMerge$TDProp=1-BuMerge$BuProp

# dev analyses
devMod<-gam(TDProp~s(interview_age,k=4)+FD+sex+RemainingTRs,data=BuMerge)
ndevMod<-gam(TDProp~FD+sex+RemainingTRs,data=BuMerge)
anovaRes<-anova.gam(devMod,ndevMod,test='Chisq')
anovaRes$`Pr(>Chi)`
```

    ## [1]           NA 1.729163e-14

``` r
# difference in adjusted r2
summary(devMod)$r.sq-summary(ndevMod)$r.sq
```

    ## [1] 0.1393841

``` r
# get residuals of model for plotting
TDProp_lm<-lm(TDProp~FD+sex+RemainingTRs,data=BuMerge)
BuMerge$TDProp_resid<-TDProp_lm$residuals+mean(BuMerge$TDProp)

ggplot(BuMerge,aes(x=(interview_age)/12,y=TDProp_resid)) + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('black'), fill = "#ed2224", alpha = .4)+ggtitle('Whole-Cortex') + geom_point(size=3,alpha=.5) + xlab("Age (Years)") +ylab("Top-Down Propagations") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), plot.title=element_text(size=30), axis.text = element_text(face="bold",size=30), axis.title = element_text(size=30, face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +scale_y_continuous(labels = percent_format(accuracy = 1))
```

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# saved at 790x810

# aaaand combat chunk on it
library(devtools)
```

    ## Loading required package: usethis

``` r
#install_github("jfortin1/neuroCombatData")
#install_github("jfortin1/neuroCombat_Rpackage")

library(neuroCombat)

# get site info
BuMerge$d_site<-droplevels(as.factor(BuMerge$site))
batch <- BuMerge$d_site
BuMerge.toharmonize <- BuMerge[,c("TDProp","FD")]
mod <- model.matrix(~BuMerge$interview_age+BuMerge$sex)
dat <- t(BuMerge.toharmonize)
hcpd.data.combat <- neuroCombat(dat=dat,mod=mod,batch=batch,eb=FALSE)
```

    ## [neuroCombat] Performing ComBat without empirical Bayes (L/S model)
    ## [neuroCombat] Found 4 batches
    ## [neuroCombat] Adjusting for  2  covariate(s) or covariate level(s)
    ## [neuroCombat] Standardizing Data across features
    ## [neuroCombat] Fitting L/S model and finding priors
    ## [neuroCombat] Adjusting the Data

``` r
dat.harmonized<-data.frame(t(hcpd.data.combat$dat.combat))

BuMerge$harmTDProp<-dat.harmonized$TDProp
BuMerge$harmFD<-dat.harmonized$FD

# get delta adjusted r^2
adjR2=summary(gam(harmTDProp~s(interview_age,k=4)+harmFD+sex,data=BuMerge))$r.sq
adjR2_noage=summary(gam(harmTDProp~harmFD+sex,data=BuMerge))$r.sq
adjR2-adjR2_noage
```

    ## [1] 0.1202557

``` r
# get anova-p
ageMod=gam(harmTDProp~s(interview_age,k=4)+harmFD+sex,data=BuMerge)
noAgeMod=gam(harmTDProp~harmFD+sex,data=BuMerge)
anovaRes=anova.gam(ageMod,noAgeMod,test='Chisq')
unlist(anovaRes$`Pr(>Chi)`)[2]
```

    ## [1] 1.987538e-12

``` r
# get residuals of model for plotting
TDProp_lm<-lm(harmTDProp~harmFD+sex,data=BuMerge)
BuMerge$TDProp_resid<-TDProp_lm$residuals+mean(BuMerge$TDProp)

ggplot(BuMerge,aes(x=(interview_age)/12,y=TDProp_resid)) + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('black'), fill = "#ed2224", alpha = .4) + geom_point(size=3,alpha=.5) + xlab("Age") +ylab("Top-Down Propagations:\n ComBat Harmonized") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=30), axis.title = element_text(size=30, face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 47, b = 0, l = 0),vjust=-4)) +scale_y_continuous(labels = percent_format(accuracy = 1))
```

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
#### show DMN Segreg

# dev analyses controlling for DMN seg
devMod_DMSC<-gam(TDProp~s(interview_age,k=4)+FD+sex+RemainingTRs+DMNsegGro,data=BuMerge)
ndevMod_DMSC<-gam(TDProp~FD+sex+RemainingTRs+DMNsegGro,data=BuMerge)

# get delta adjusted r^2
adjR2=summary(gam(TDProp~s(interview_age,k=4)+FD+sex+RemainingTRs+DMNsegGro,data=BuMerge))$r.sq
adjR2_noage=summary(gam(TDProp~FD+sex+RemainingTRs+DMNsegGro,data=BuMerge))$r.sq
adjR2-adjR2_noage
```

    ## [1] 0.1414908

``` r
# get p-value
anovaRes<-anova.gam(devMod_DMSC,ndevMod_DMSC,test='Chisq')
anovaRes$`Pr(>Chi)`
```

    ## [1]           NA 1.201915e-14

``` r
# DMN model - group
summary(gam(DMNsegGro~s(interview_age,k=4)+FD+sex+RemainingTRs+TDProp,data=BuMerge))
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## DMNsegGro ~ s(interview_age, k = 4) + FD + sex + RemainingTRs + 
    ##     TDProp
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  -5.559e-02  3.463e-02  -1.605   0.1093  
    ## FD           -2.733e-02  6.358e-02  -0.430   0.6675  
    ## sexM          6.061e-04  9.550e-04   0.635   0.5260  
    ## RemainingTRs -4.970e-06  2.009e-06  -2.474   0.0138 *
    ## TDProp        7.144e-02  6.799e-02   1.051   0.2941  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                    edf Ref.df     F  p-value    
    ## s(interview_age) 1.708  2.074 9.684 6.78e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.102   Deviance explained = 11.6%
    ## GCV = 8.8834e-05  Scale est. = 8.7298e-05  n = 388

``` r
devMod_DMgro<-gam(DMNsegGro~s(interview_age,k=4)+FD+sex+RemainingTRs+TDProp,data=BuMerge)
ndevMod_DMgro<-gam(DMNsegGro~FD+sex+RemainingTRs+TDProp,data=BuMerge)

# get delta adjusted r^2
adjR2=summary(gam(DMNsegGro~s(interview_age,k=4)+FD+sex+RemainingTRs+TDProp,data=BuMerge))$r.sq
adjR2_noage=summary(gam(DMNsegGro~FD+sex+RemainingTRs+TDProp,data=BuMerge))$r.sq
adjR2-adjR2_noage
```

    ## [1] 0.04510866

``` r
# get p-value
anovaRes<-anova.gam(devMod_DMgro,ndevMod_DMgro,test='Chisq')
anovaRes$`Pr(>Chi)`
```

    ## [1]           NA 3.142891e-05

``` r
# and plotting dmn seg directly
DMN_lm<-lm(DMNseg~FD+sex+RemainingTRs+TDProp,data=BuMerge)
BuMerge$DMN_lm_resid<-DMN_lm$residuals+mean(BuMerge$DMNseg)

ggplot(BuMerge,aes(x=(interview_age)/12,y=DMN_lm_resid)) + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('black'), fill = "#ed2224", alpha = .4,size=3) + geom_point(size=5,alpha=.5) + xlab("Age") +ylab("DMN External \n Connectivity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=30), axis.title = element_text(size=30, face="bold"), axis.title.y = element_text(margin = margin(t = 0, r = 47, b = 0, l = 0),vjust=-4))
```

![](Group_level_analysis_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
# only written once, no need to overwrite
#write.table(masterdf$SubjID,'~/G600_TRs.txt',quote = F,col.names = F,row.names = F)
```

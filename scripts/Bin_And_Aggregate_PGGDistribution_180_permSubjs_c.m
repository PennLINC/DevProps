% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
Subjs=readtable('~/PWs/rs_cmatch_subs.csv')
%initialize output array
% 18 bins for angular distances because 0 and 180 are endpoints (19 numbers)
OutDf=zeros(height(Subjs),18);
OutDf_c=zeros(height(Subjs),18);
% 17 bins for modes because there is no 0th bin - 1 and 18 are endpoints (18 numbers)
OutDf_modes=zeros(1,17);
% edges to apply in discretize
edges = 0:10:180;
modeedges = 1:18;
%%% prepare mask
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
% set diff to get those needed
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_AngDistMat4.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_AngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% discretize : get all TR pairs from non mw-indices
		Disc_FAngles_R=histcounts(AngsR(:,g_noMW_combined_R),edges);
		Disc_FAngles_L=histcounts(AngsL(:,g_noMW_combined_L),edges);
		OutDf(s,:)=Disc_FAngles_L;
		OutDf(s,:)=Disc_FAngles_R;
	end
	% same for _carit, but drop the conditional
        % load in subj's distr
        Angs=load([fp '/' Subj{:} '_AngDistMat4_c.mat']);
        AngsL=Angs.AngDist.gLeft;
        AngsR=Angs.AngDist.gRight;
        % discretize : get all TR pairs from non mw-indices
        Disc_FAngles_R=histcounts(AngsR(:,g_noMW_combined_R),edges);
        Disc_FAngles_L=histcounts(AngsL(:,g_noMW_combined_L),edges);
        OutDf_c(s,:)=Disc_FAngles_L;
        OutDf_c(s,:)=Disc_FAngles_R;
end

% stack both of them to randomly sample from
OutDf_master=vertcat(OutDf,OutDf_c);

% initialize difference mat
difMat=zeros(1000,18);
% old and young subj groups are 127 and 132 subjs, respectively. Get random groups of equal size over perms
% now we permute 1000 times
for p=1:1000
	% permute subjseq to get ordering to drawfrom
	subjseqP=randperm(562);
	% get pseudo old group
	pseuRest=subjseqP(1:281);
	% get pseudo young group
	pseuTask=subjseqP(282:562);
	% index out
	pseuRestmat=OutDf_master(pseuRest,:);
	pseuTaskmat=OutDf_master(pseuTask,:);
	% normalize to percentages
	totalR=sum(pseuRestmat,2);
	normR=pseuRestmat./totalR;
	totalT=sum(pseuTaskmat,2);
	normT=pseuTaskmat./totalT;
	% get mean
	mR=mean(normR);
	mT=mean(normT);
	% difference in percentages
	difMat(p,:)=mT-mR;
	p
end	

%%% get observed difference 
% normalize and tag onto end of df
totalR=sum(OutDf,2);
totalT=sum(OutDf_c,2);
normR=OutDf./totalR;
normT=OutDf_c./totalT;
% get mean
mR=mean(normR);
mT=mean(normT);

% tag onto very end of difMat
difMat(1001,:)=mT-mR;

% write out as table
fn=['/cbica/projects/pinesParcels/results/PWs/Task_v_rest_AngDistDif_permuted.csv'];
writetable(table(difMat),fn)

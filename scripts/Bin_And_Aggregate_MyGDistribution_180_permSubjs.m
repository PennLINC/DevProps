% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
%%%% Age Tertiles as well
Subjs=readtable('~/PWs/rs_subs.csv')
%initialize output array
% 18 bins for angular distances because 0 and 180 are endpoints (19 numbers)
OutDf=zeros(height(Subjs),18);
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
gLPGfp=['~/data/hcpd.mm_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['~/data/hcpd.mm_R_3k.func.gii'];
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
	if isfile([fp '/' Subj{:} '_Myelin_AngDistMat4.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_Myelin_AngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% discretize : get all TR pairs from non mw-indices
		Disc_FAngles_R=histcounts(AngsR(:,g_noMW_combined_R),edges);
		Disc_FAngles_L=histcounts(AngsL(:,g_noMW_combined_L),edges);
		OutDf(s,:)=Disc_FAngles_L;
		OutDf(s,:)=Disc_FAngles_R;
	end
end

% initialize difference mat
difMat=zeros(1000,18);
% old and young subj groups are 127 and 132 subjs, respectively. Get random groups of equal size over perms
% now we permute 1000 times
for p=1:1000
	% permute subjseq to get ordering to drawfrom
	subjseqP=randperm(388);
	% get pseudo old group
	pseuOld=subjseqP(1:132);
	% get pseudo young group
	pseuYoung=subjseqP(133:259);
	% index out
	pseuOldmat=OutDf(pseuOld,:);
	pseuYoungmat=OutDf(pseuYoung,:);
	% normalize to percentages
	totalOl=sum(pseuOldmat,2);
	normOl=pseuOldmat./totalOl;
	totalYo=sum(pseuYoungmat,2);
	normYo=pseuYoungmat./totalYo;
	% get mean
	mOl=mean(normOl);
	mYo=mean(normYo);
	% difference in percentages
	difMat(p,:)=mOl-mYo;
	p
end	

% get observed difference 
youngs=readtable('~/PWs/young_subs.csv');
olds=readtable('~/PWs/old_subs.csv');

% convert Subjs to table specifically for indexing
T=table(Subjs,'RowNames',Subjs{:,2});
TrueYoung=T(youngs{:,2},:);
TrueOld=T(olds{:,2},:);

% ensure we are pulling the correct subj IDs
if height(unique(vertcat(youngs(:,2),olds(:,2)))) ~= length(pseuOld)+length(pseuYoung)
	message='Young and Old observed vs permuted lengths do not match'
	throw(message)
end

% pull out indices
TrueOld=TrueOld{:,1};
IndicesOld=TrueOld{:,1};
TrueYoung=TrueYoung{:,1};
IndicesYoung=TrueYoung{:,1};

% pull out values
TrueOldmat=OutDf(str2double(IndicesOld),:);
TrueYoungmat=OutDf(str2double(IndicesYoung),:);

% normalize and tag onto end of df
totalOl=sum(TrueOldmat,2);
totalYo=sum(TrueYoungmat,2);
normOl=TrueOldmat./totalOl;
normYo=TrueYoungmat./totalYo;
% get mean
mOl=mean(normOl);
mYo=mean(normYo);

% tag onto very end of difMat
difMat(1001,:)=mOl-mYo;

% write out as table
fn=['/cbica/projects/pinesParcels/results/PWs/Old_v_young_AngDistDif_MyG_permuted.csv'];
writetable(table(difMat),fn)

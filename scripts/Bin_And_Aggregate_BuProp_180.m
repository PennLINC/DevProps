% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
%%%%%%% MODULAR: WILL NEED RS, RS MATCHED FROM _C, and _C SUBJ LIST
%%%% Age Tertiles as well
Subjs=readtable('~/PWs/rs_subs.csv')
%initialize output array
% 18 bins for angular distances because 0 and 180 are endpoints (19 numbers)
OutDf=zeros(1,height(Subjs));
% 17 bins for modes because there is no 0th bin - 1 and 18 are endpoints (18 numbers)
% edges to apply in discretize
edges = 0:10:180;
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
		% initialize facewise vectors for 1-18 mode for 
		facePropL=zeros(1,length(g_noMW_combined_L));
		facePropR=zeros(1,length(g_noMW_combined_R));
		% for each face: get mode over TRs
		for l=1:length(g_noMW_combined_L);
			binnedFaceDirs=histcounts(AngsL(:,g_noMW_combined_L(l)),edges);
			BUs=sum(binnedFaceDirs(1:9));
			TDs=sum(binnedFaceDirs(10:18));
			facePropL(l)=BUs/(BUs+TDs);	
		end
		% for right
		for r=1:length(g_noMW_combined_R);
                        binnedFaceDirs=histcounts(AngsR(:,g_noMW_combined_R(r)),edges);
                        BUs=sum(binnedFaceDirs(1:9));
                        TDs=sum(binnedFaceDirs(10:18));
                        facePropR(r)=BUs/(BUs+TDs);
		end
		% count of prop BU for this subj
		OutDf(s)=mean([facePropL facePropR]);
	end
end

% maybe write as csv: as table ya know
fn=['/cbica/projects/pinesParcels/results/PWs/rs_subsBuProp.csv'];
writetable(table(OutDf),fn)


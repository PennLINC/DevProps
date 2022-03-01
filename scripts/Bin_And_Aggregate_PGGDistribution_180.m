% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
%%%%%%% MODULAR: WILL NEED RS, RS MATCHED FROM _C, and _C SUBJ LIST
%%%% Perhaps Age Tertiles as well
Subjs=readtable('~/PWs/old_subs.csv')
%initialize output array
OutDf=zeros(1,18);
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
	if isfile([fp '/' Subj{:} '_AngDistMat4_c.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_AngDistMat4_c.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% discretize : get all TR pairs from non mw-indices
		Disc_FAngles_R=histcounts(AngsR(:,g_noMW_combined_R),edges);
		Disc_FAngles_L=histcounts(AngsL(:,g_noMW_combined_L),edges);
		OutDf=OutDf+Disc_FAngles_L;
		OutDf=OutDf+Disc_FAngles_R;	
	end
end

% maybe write as csv: as table ya know
fn=['/cbica/projects/pinesParcels/results/PWs/old_subs_AngDistHist.csv'];
writetable(table(OutDf),fn)



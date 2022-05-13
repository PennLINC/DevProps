% add paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in subj list
Subjs=readtable('~/PWs/rs_subs.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mw mask
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4';
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
% use native freesurfer command for mw mask indices
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% init df
outDF=zeros(height(Subjs),(length(g_noMW_combined_L)+length(g_noMW_combined_R)));
% for each subj
for s=1:height(Subjs)
	s
	Subj=table2array(Subjs(s,2));
	subj=Subj{:};
	% load in
	fpL=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
	dataL=load(fpL);
	fpR=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
	dataR=load(fpR);
	% dist from PGG in col 8, valid faces in rows specified (non-mw)
	PGGd_L=cell2mat(dataL.OutDf_L(g_noMW_combined_L,8));
	PGGd_R=cell2mat(dataR.OutDf_R(g_noMW_combined_R,8));
	% convert to degrees
	PGGd_L=rad2deg(PGGd_L);
	PGGd_R=rad2deg(PGGd_R);
	% convert to 0-360 range
	PGGd_L=abs(PGGd_L);
	PGGd_R=abs(PGGd_R);
	% convert to 0-180 range
	% we will subtract pi from every value so as to center the distribution on pi for reflecting
	PGGd_L_min=PGGd_L-180;
	PGGd_R_min=PGGd_R-180;
	% then we reflect it with abs()
	PGGd_L_minabs=abs(PGGd_L_min);
	PGGd_R_minabs=abs(PGGd_R_min);
	% we've successfully collapsed half of our distribution, but now our distribution is backwards: lowest values are high
	PGGd_L_minabsneg=PGGd_L_minabs*-1;
	PGGd_R_minabsneg=PGGd_R_minabs*-1;
	%  now our distribution is collapsed and in the right order, but on the wrong side of the tracks. restore the initial pi we subtracted.
	PGGd_L_180=PGGd_L_minabsneg+180;
	PGGd_R_180=PGGd_R_minabsneg+180;
	% add to df
	outDF(s,:)=vertcat(PGGd_L_180,PGGd_R_180);
end
% saveout in r-friendly format

% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
%%%%%%% MODULAR: WILL NEED RS, RS MATCHED FROM _C, and _C SUBJ LIST
%%%% Age Tertiles as well
Subjs=readtable('~/PWs/rs_subs.csv')
%initialize output array
% 18 bins for angular distances because 0 and 180 are endpoints (19 numbers)
OutDf=zeros(1,18);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% and these are for group-level modes at each face
OutDf_face_modesL=zeros(height(Subjs),length(g_noMW_combined_L));
OutDf_face_modesR=zeros(height(Subjs),length(g_noMW_combined_R));
% add partner DFs for mode relative prominence
OutDf_face_modesL_Prom=zeros(height(Subjs),length(g_noMW_combined_L));
OutDf_face_modesR_Prom=zeros(height(Subjs),length(g_noMW_combined_R));

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_curvAngDistMat4.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_curvAngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% discretize : get all TR pairs from non mw-indices
		Disc_FAngles_R=histcounts(AngsR(:,g_noMW_combined_R),edges);
		Disc_FAngles_L=histcounts(AngsL(:,g_noMW_combined_L),edges);
		OutDf=OutDf+Disc_FAngles_L;
		OutDf=OutDf+Disc_FAngles_R;
		% initialize facewise vectors for 1-18 mode for 
		faceModesL=zeros(1,length(g_noMW_combined_L));
		faceModesR=zeros(1,length(g_noMW_combined_R));
		faceModesL_Prom=zeros(1,length(g_noMW_combined_L));
		faceModesR_Prom=zeros(1,length(g_noMW_combined_R));
		% for each face: get mode over TRs
		for l=1:length(g_noMW_combined_L);
			binnedFaceDirs=histcounts(AngsL(:,g_noMW_combined_L(l)),edges);
			[M,I]=max(binnedFaceDirs);
			faceModesL(l)=I;
			%%% find prominence from non-adjacent modes
			% make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
			tmpBFDvec=zeros(1,22);
			tmpBFDvec(3:20)=binnedFaceDirs;
			% I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
			tmpBFDvec((I):(I+4))=[];
			% get updated binnedFaceDirs: those non-adjacent
			[M2,I2]=max(tmpBFDvec);
			% get relative prominence of secondary peak
			faceModesL_Prom(l)=1-(M2/M);
		end
		% for right
		for r=1:length(g_noMW_combined_R);
                        binnedFaceDirs=histcounts(AngsR(:,g_noMW_combined_R(r)),edges);
                        [M,I]=max(binnedFaceDirs);
                        faceModesR(r)=I;
			%%% find prominence from non-adjacent modes
                        % make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
                        tmpBFDvec=zeros(1,22);
                        tmpBFDvec(3:20)=binnedFaceDirs;
                        % I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
                        tmpBFDvec((I):(I+4))=[];
                        % get updated binnedFaceDirs: those non-adjacent
                        [M2,I2]=max(tmpBFDvec);
                        % get relative prominence of secondary peak
                        faceModesR_Prom(r)=1-(M2/M);
                end
		% get histcounts of vector of modes for each face (collapsed over TRs)	
		Disc_modes_FAngles_L=histcounts(faceModesL,modeedges);
		Disc_modes_FAngles_R=histcounts(faceModesR,modeedges);
		% stack in to group-level
		OutDf_modes=OutDf_modes+Disc_modes_FAngles_L;
                OutDf_modes=OutDf_modes+Disc_modes_FAngles_R;
		% get mode AT each face, to be once again moded for plotting group-level modes
		OutDf_face_modesL(s,:)=faceModesL;
		OutDf_face_modesR(s,:)=faceModesR;
		% save mode prominence for this subj
		OutDf_face_modesL_Prom(s,:)=faceModesL_Prom;
		OutDf_face_modesR_Prom(s,:)=faceModesR_Prom;	
	end
end

% take mode of modes for faces
OutDf_face_modeL=mode(OutDf_face_modesL,1);
OutDf_face_modeR=mode(OutDf_face_modesR,1);
% get mean prominence
OutDf_face_modeL_prom=mean(OutDf_face_modesL_Prom,1);
OutDf_face_modeR_prom=mean(OutDf_face_modesR_Prom,1);

fn=['/cbica/projects/pinesParcels/results/PWs/rs_subs_facewiseMode_curvL.csv'];
writetable(table(OutDf_face_modeL,OutDf_face_modeL_prom),fn)
fn=['/cbica/projects/pinesParcels/results/PWs/rs_subs_facewiseMode_curvR.csv'];
writetable(table(OutDf_face_modeR,OutDf_face_modeR_prom),fn)


% maybe write as csv: as table ya know
fn=['/cbica/projects/pinesParcels/results/PWs/rs_subs_AngDistHist_curv.csv'];
writetable(table(OutDf),fn)

% write out modes sep.
fn=['/cbica/projects/pinesParcels/results/PWs/rs_subs_AngDistHist_modedOverTRs_curv.csv'];
writetable(table(OutDf_modes),fn)

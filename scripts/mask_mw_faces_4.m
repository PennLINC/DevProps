function mask_mw_faces(subj)
% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in angular distances
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat4.mat'];
data=load(AngDistFP)

%%% Load in surface data
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

% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);

% all faces that touch a medial wall vertex to be masked
%MW_f1_L=find(ismember(F_L(:,1),mw_v_l));
%MW_f2_L=find(ismember(F_L(:,2),mw_v_l));
%MW_f3_L=find(ismember(F_L(:,3),mw_v_l));
% rh
%MW_f1_R=find(ismember(F_R(:,1),mw_v_r));
%MW_f2_R=find(ismember(F_R(:,2),mw_v_r));
%MW_f3_R=find(ismember(F_R(:,3),mw_v_r));
% inclusive - mask if face involves ANY mw vertices
%MW_combined_L=union(MW_f1_L,MW_f2_L);
%MW_combined_L=union(MW_combined_L,MW_f3_L);
% now for right hemisphere
%MW_combined_R=union(MW_f1_R,MW_f2_R);
%MW_combined_R=union(MW_combined_R,MW_f3_R);
% to be combined with mask vector derived below from PGGs


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

% mask angular distances
gAngD_L=data.AngDist.gLeft;
gAngD_R=data.AngDist.gRight;
gAngD_L_masked=gAngD_L(:,g_noMW_combined_L);
gAngD_R_masked=gAngD_R(:,g_noMW_combined_R);

%%% mask PG
% convert PG to faces, average value across three vertices defining face
% for f in faces
% vert 1 = faces(f,1)
% vert 2 = faces(f,2)
% vert 3 = faces(f,3)
% pg val = mean(PG(vert1),PG(vert2),PG(vert3)
%end


% write in R-friendly format
writetable(table(gAngD_L_masked),['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_gAngDist_Masked4_L.csv']);
writetable(table(gAngD_R_masked),['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_gAngDist_Masked4_R.csv']);


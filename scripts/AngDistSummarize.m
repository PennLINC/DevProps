function [numBU, numTD] = AngDistSummarize(subj)
% return number of frames with mean Ang Dist > 120, # < 60

%%%% Load in angular distances
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat4.mat'];
data=load(AngDistFP);

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
MW_f1_L=find(ismember(F_L(:,1),mw_v_l));
MW_f2_L=find(ismember(F_L(:,2),mw_v_l));
MW_f3_L=find(ismember(F_L(:,3),mw_v_l));
% rh
MW_f1_R=find(ismember(F_R(:,1),mw_v_r));
MW_f2_R=find(ismember(F_R(:,2),mw_v_r));
MW_f3_R=find(ismember(F_R(:,3),mw_v_r));
% inclusive - mask if face involves ANY mw vertices
MW_combined_L=union(MW_f1_L,MW_f2_L);
MW_combined_L=union(MW_combined_L,MW_f3_L);
% now for right hemisphere
MW_combined_R=union(MW_f1_R,MW_f2_R);
MW_combined_R=union(MW_combined_R,MW_f3_R);
% to be combined with mask vector derived below from PGGs

% MASK WHERE PGG = 0: individ AND group
% load in GROUP PG
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

% get inverse for indexing : faces that ARE NOT touching mW verts
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

% mask angular distances
gAngD_L=data.AngDist.gLeft;
gAngD_R=data.AngDist.gRight;
gAngD_L_masked=gAngD_L(:,g_noMW_combined_L);
gAngD_R_masked=gAngD_R(:,g_noMW_combined_R);

% get mean angular distance per TR
lmAd=mean(gAngD_L_masked,2);
rmAd=mean(gAngD_R_masked,2);

numBU=sum(rmAd<60);
numTD=sum(rmAd>120);


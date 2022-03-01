%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Save out the distribution of vector magnitudes (PGG) for mask-masking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

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
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

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

% extract face-wise vector cartesian vector components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);

% sum vectors - absolute values
vsumL=sum(abs(gPGg_L),2);
vsumR=sum(abs(gPGg_R),2);

% hypotenuse of vectors
vHL=sqrt((gPGx_L.^2)+(gPGy_L.^2)+(gPGz_L.^2));
vHR=sqrt((gPGx_R.^2)+(gPGy_R.^2)+(gPGz_R.^2));

% generate mw mask
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

% apply mw mask
vsumL=vsumL(g_noMW_combined_L);
vsumR=vsumR(g_noMW_combined_R);
vHL=vHL(g_noMW_combined_L);
vHR=vHR(g_noMW_combined_R);

% save for R
writetable(table(vsumL),'~/data/vsumL.csv')
writetable(table(vsumR),'~/data/vsumR.csv')
writetable(table(vHL),'~/data/vHL.csv')
writetable(table(vHR),'~/data/vHR.csv')


% read in spun PGs, convert to single cifti for more rapid null-direction testing (1 file load rather than 1000 in subj runs)

% needed paths - note cifti matlab not usual toolbox path
addpath(genpath('/cbica/projects/pinesParcels/cifti-matlab'))

% load in spins
sp_fp='/cbica/projects/pinesParcels/results/aggregated_data/PGPermuts_fs4.mat';
sp=load(sp_fp)

% load template cifti
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% replace cdata

cifti_create_dscalar_from_template

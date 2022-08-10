% load in PGs, save em out as r-friendly csvs
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in vertices PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.lh.fsaverage5.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH_V=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.rh.fsaverage5.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH_V=gRPGf.cdata(:,1);

writetable(table(gPG_LH_V),'~/results/PWs/vertPG_left.csv')
writetable(table(gPG_RH_V),'~/results/PWs/vertPG_right.csv')

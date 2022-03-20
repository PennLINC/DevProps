% make a .mat version of smoothed PG1 because compiled matlablab code can't get "gifti" to work without an .xml error
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/lh_fs4.avg_curv.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/rh_fs4.avg_curv.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/fs4curv.mat','gpg')

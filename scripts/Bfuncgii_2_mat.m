% make a .mat version of smoothed PG1 because compiled matlablab code can't get "gifti" to work without an .xml error
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/beta3k_L.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/beta3k_R.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/gB_fs4.mat','gpg')

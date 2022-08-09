% make a .mat version of myelin map because compiled matlablab code can't get "gifti" to work without an .xml error
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in GROUP PG
gLPGfp=['~/data/hcpd.mm_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['~/data/hcpd.mm_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/mypg_fs4.mat','gpg')

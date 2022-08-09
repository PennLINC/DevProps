% make a .mat version of smoothed PG1 because compiled matlablab code can't get "gifti" to work without an .xml error
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in young tertile PG
gLPGfp=['~/data/hcpd.youngPG_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['~/data/hcpd.youngPG_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/ypg_fs4.mat','gpg')

% load in mid tertile PG
gLPGfp=['~/data/hcpd.midPG_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['~/data/hcpd.midPG_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/mipg_fs4.mat','gpg')

% load in old tertile PG
gLPGfp=['~/data/hcpd.oldPG_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['~/data/hcpd.oldPG_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% saveout
gpg=struct;
gpg.gPG_RH=gPG_RH;
gpg.gPG_LH=gPG_LH;
save('~/data/opg_fs4.mat','gpg')


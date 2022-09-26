% make a .mat version of smoothed PG1 because compiled matlablab code can't get "gifti" to work without an .xml error
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in pial surface for most accurate rep. of cortical A->P
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.pial'];
surfR = [SubjectsFolder '/surf/rh.pial'];
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% initialize compile-friendly struct, same struct naming as PG
gpg=struct;
% second coordinate set in vx_* is a->p, tested with vis_vertvec
gpg.gPG_LH=vx_l(:,2);
gpg.gPG_RH=vx_r(:,2);
% saveout
save('~/data/fs4_ap.mat','gpg')

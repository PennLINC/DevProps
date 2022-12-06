% addpaths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
WorkingFolder = '/cbica/projects/pinesParcels/data/SingleParcellation/SingleAtlas_Analysis';
% set outdir
outdir='/cbica/projects/pinesParcels/results/aggregated_data/';

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
pg1=[gpg.gPG_LH' gpg.gPG_RH'];

% set output file name
outFn=strcat('/cbica/projects/pinesParcels/results/aggregated_data/PAPermuts_fs4.mat');

% load in mask (SNR Mask)
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
grad_lh(mwIndVec_l)=100;
grad_rh(mwIndVec_r)=100;

% write them out as a csv for spin test to deal with		
writetable(table(grad_lh),[outdir 'PAfs4_L.csv'],'WriteVariableNames',0);
writetable(table(grad_rh),[outdir 'PAfs4_R.csv'],'WriteVariableNames',0);
% create permutations, save out to outFn
SpinPermuFS([outdir 'PAfs4_L.csv'], [outdir 'PAfs4_R.csv'], 1000, outFn);

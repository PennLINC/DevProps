% addpaths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
WorkingFolder = '/cbica/projects/pinesParcels/data/SingleParcellation/SingleAtlas_Analysis';
% set outdir
outdir='/cbica/projects/pinesParcels/results/aggregated_data/';

% get gradient map
pgl = gifti(['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii']);
pgr = gifti(['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii']);
grad_lh = pgl.cdata(:,1);
grad_rh = pgr.cdata(:,1);
pg1=[grad_lh' grad_rh'];

% set output file name
outFn=strcat('/cbica/projects/pinesParcels/results/aggregated_data/PGPermuts_fs4.mat');

% load in mask (SNR Mask)
surfML = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
grad_lh(mwIndVec_l)=100;
grad_rh(mwIndVec_r)=100;

% write them out as a csv for spin test to deal with		
writetable(table(grad_lh),[outdir 'gradfs4_L.csv'],'WriteVariableNames',0);
writetable(table(grad_rh),[outdir 'gradfs4_R.csv'],'WriteVariableNames',0);
% create permutations, save out to outFn
SpinPermuFS([outdir 'gradfs4_L.csv'], [outdir 'gradfs4_R.csv'], 1000, outFn);

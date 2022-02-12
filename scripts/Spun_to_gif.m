

% needed paths - slightly different than standard toolbox fp, had to follow an obscure readthedocs
addpath(genpath('/cbica/projects/pinesParcels/gifti'))

% load in spins
sp_fp='/cbica/projects/pinesParcels/results/aggregated_data/PGPermuts_fs4.mat';
sp=load(sp_fp);

% load in template file
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
t_lh=gifti(gLPGfp);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
t_rh=gifti(gRPGfp);

% extract spun surface
data_lh=sp.bigrotl;
data_rh=sp.bigrotr;

% template from https://enigma-toolbox.readthedocs.io/en/latest/pages/18.export/index.html
data = gifti(CT_schaefer_200_c69);
data_lh = data; data_lh.cdata = data_lh.cdata(1:end/2);
data_rh = data; data_rh.cdata = data_rh.cdata(end/2+1:end);

% Define output path and filenames
dpath = what('import_export'); dpath = dpath.path;
fname_lh = 'lh.schaefer_200_c69_thickness.gii'
fname_rh = 'rh.schaefer_200_c69_thickness.gii'

% Export data as GIfTI / .gii
savegifti(data_lh, [dpath, fname_lh], 'Base64Binary');
savegifti(data_rh, [dpath, fname_rh], 'Base64Binary');

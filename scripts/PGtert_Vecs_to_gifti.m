% add paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in template data with gifti
template_l=gifti('/cbica/projects/pinesParcels/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii');
template_r=gifti('/cbica/projects/pinesParcels/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii');

% load in pg values
ypg=readtable('/cbica/projects/pinesParcels/results/PWs/FC/young_PG1.csv','ReadVariableNames',false);
mipg=readtable('/cbica/projects/pinesParcels/results/PWs/FC/mid_PG1.csv','ReadVariableNames',false);
opg=readtable('/cbica/projects/pinesParcels/results/PWs/FC/old_PG1.csv','ReadVariableNames',false);


% run over young, old and mid


% replace template data with pg values
template_l.cdata=ypg{1:10242,1};
template_r.cdata=ypg{10243:20484,1};

% saveout
save(template_l,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_young_L.func.gii')
save(template_r,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_young_R.func.gii')

% replace template data with pg values
template_l.cdata=mipg{1:10242,1};
template_r.cdata=mipg{10243:20484,1};

% saveout
save(template_l,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_mid_L.func.gii')
save(template_r,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_mid_R.func.gii')

% replace template data with pg values
template_l.cdata=opg{1:10242,1};
template_r.cdata=opg{10243:20484,1};

% saveout
save(template_l,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_old_L.func.gii')
save(template_r,'/cbica/projects/pinesParcels/data/princ_gradients/hcpd_old_R.func.gii')


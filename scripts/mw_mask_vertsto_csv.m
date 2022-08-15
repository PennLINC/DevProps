% add paths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% read in mwvecs
surfML = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);

% binary mask
mw_L=ones(1,10242);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,10242);
mw_R(mwIndVec_r)=0;

%saveout 
writetable(table(mw_L),'~/data/mw_boolean_l.csv')
writetable(table(mw_R),'~/data/mw_boolean_r.csv')

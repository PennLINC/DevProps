% add paths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% read in mwvecs
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);

% binary mask
mw_L=ones(1,2562);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,2562);
mw_R(mwIndVec_r)=0;

% check for .2% of missed vertices from ->fs4 transform


% tertile list
tertiles={'young','old','mid'}
% load in medial walls
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "isn't medial wall" vector for vertices
mw_L=ones(1,2562);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,2562);
mw_R(mwIndVec_r)=0;
% for each tertile
tertile=tertiles{1}
subjs=fileread(['~/PWs/' tertile '_subs.txt'])
% matlab is picky eater
subjs=strsplit(subjs);
subjs=subjs(1:(length(subjs)-1));
s=1
subj=subjs(s)
vw_ts_l_p=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_3k.mgh'];
vw_ts_r_p=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_3k.mgh'];
% this is for sure some fussy matlab bs
vw_ts_l_p=strcat(vw_ts_l_p{1},vw_ts_l_p{2},vw_ts_l_p{3},vw_ts_l_p{4},vw_ts_l_p{5});
vw_ts_r_p=strcat(vw_ts_r_p{1},vw_ts_r_p{2},vw_ts_r_p{3},vw_ts_r_p{4},vw_ts_r_p{5});
vw_ts_l=MRIread(vw_ts_l_p);
vw_ts_r=MRIread(vw_ts_r_p);
% series length
vw_ts_l_SIZE=size(vw_ts_l.vol);
vw_ts_l_LEN=vw_ts_l_SIZE(4)
vw_ts_l=vw_ts_l.vol;
vw_ts_r=vw_ts_r.vol;
% stacking matrices so vertex number is doubled (not timepoints)
vw_ts_both=[vw_ts_l vw_ts_r];
vw_ts_bothrw=zeros(vw_ts_l_LEN,5124);
% calculate FC
for x=1:length(vw_ts_bothrw)
    vw_ts_bothrw(:,x)=vw_ts_both(1,x,1,:);
end
% bigass connectivity matrix, takes 5 seconds or so to calc
ba_conmat=corrcoef(vw_ts_bothrw);
% left
ba_conmat_L=ba_conmat(1:2562,1:2562);
% right
ba_conmat_R=ba_conmat(2563:5124,2563:5124);
% use mask
ba_conmat_L_m=ba_conmat_L(logical(mw_L),logical(mw_L));
ba_conmat_R_m=ba_conmat_R(logical(mw_R),logical(mw_R));
% get remaining .2% of vertices that are NaNs
NansL=isnan(ba_conmat_L_m(1,:));
NansR=isnan(ba_conmat_R_m(1,:));
% add nanspots to mask
mw_L(logical(mw_L))=mw_L(logical(mw_L))-NansL;
mw_R(logical(mw_R))=mw_R(logical(mw_R))-NansR;

%saveout 
writetable(table(mw_L),'~/data/mw_boolean_lfs4.csv')
writetable(table(mw_R),'~/data/mw_boolean_rfs4.csv')

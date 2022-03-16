function us = OpFl_Sph(subj)


% f should be data in {----}, both in the same order as vx_l (and vx_r) below. 
% adapted from Decomposition of Optical Flow on the Sphere, by Kirisits, Lang, and Scherzer (2014)
% https://www.csc.univie.ac.at/paper/KirLanSch14.pdf

% set to run independently on each pair of temporally adjacent frames... yikes

%%% squeeze(hemiTS.mgh.vol) before input!

%%%%%%%%%%%%%%%%%%%% Set parameters.
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles 
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
% load in TRs_l and TRs_r
TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_10k.mgh'];
TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_10k.mgh'];
% filepaths to files
TRs_lf=MRIread(TRs_lfp);
TRs_rf=MRIread(TRs_rfp);
% files to data
TRs_l=squeeze(TRs_lf.vol);
TRs_r=squeeze(TRs_rf.vol);
% load in fsaverage5 faces and vertices
%addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage5';
% for surface data
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topology
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l)
sizeInDr=size(TRs_r)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
	fl.TRs{TRP}=TRs_l(:,TRP);
end

disp('converting right hemi to struct')

% right hemi
fr=struct;
% populate struct
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r(:,TRP);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow...');

% loop over each TR-Pair: 1 fewer pair than number of TRs

for TRP=1:(TR_n-1);
% Compute decomposition.
tic;
% pull out adjacent frames
u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
toc;
% throw u into struct
us.vf_left{TRP}=u;
% now right hemi
u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
toc;
% throw u into struct
us.vf_right{TRP}=u;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs5.mat'],'us')

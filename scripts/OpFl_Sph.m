function u = OpFl_Sph(f)
% f should be data in two frames, {1} and {2}, both in the same order as vx_l (and vx_r) below. 
% adapted from Decomposition of Optical Flow on the Sphere, by Kirisits, Lang, and Scherzer (2014)
% https://www.csc.univie.ac.at/paper/KirLanSch14.pdf

% set to run independently on each pair of temporally adjacent frames... yikes

%%% squeeze(hemiTS.mgh.vol) before input!

% Set parameters.
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles 
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative

% load in fsaverage4 faces and vertices
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4';
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
vx_l(numV+1:end, :) = Normalize(vx_l(numV+1:end, :));

% right
numV=length(vx_r);
vx_l(numV+1:end, :) = Normalize(vx_l(numV+1:end, :));

% Compute decomposition.
disp('Computing optical flow...');
tic;
u = of(N, faces_l, vx_l, f{1}, f{2}, h, alpha, s);
toc;

% consider resample to fsavergage 5 after frame-by-frame aggregations are done



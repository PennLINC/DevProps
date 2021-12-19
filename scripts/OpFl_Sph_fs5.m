function us = OpFl_Sph_fs5(TRs_l,TRs_r)


% f should be data in {----}, both in the same order as vx_l (and vx_r) below. 
% adapted from Decomposition of Optical Flow on the Sphere, by Kirisits, Lang, and Scherzer (2014)
% https://www.csc.univie.ac.at/paper/KirLanSch14.pdf

% set to run independently on each pair of temporally adjacent frames... yikes

%%% convert to mgh before input!
%%% squeeze(hemiTS.mgh.vol) before input!

%%%%%%%%%%%%%%%%%%%% Set parameters.
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles 
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
% load in fsaverage4 faces and vertices
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5';
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

%disp('converting right hemi to struct')

% right hemi
%fr=struct;
% populate struct
%for TRP=1:TR_n;
%        fr.TRs{TRP}=TRs_r(:,TRP);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow...');

% loop over each TR-Pair: 1 fewer pair than number of TRs

for TRP=1:(TR_n-1);
TRP
% Compute decomposition.
tic;
% pull out adjacent frames
u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
toc;
% throw u into struct
us.vf_left{TRP}=u;
% now right hemi
%u = of(N, faces_l, vx_l, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
%toc;
% throw u into struct
%us.vf_right{TRP}=u;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('~/dropbox/OpFl_Sph_fs5_test.mat','us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualize
% extract medial wall mask
% LHmw=read_label([],'/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/lh.Medial_wall');
% LHmw_Inds=LHmw(:,1);
% might be a +1 in indexing begins at 0 again
% LHmw_Inds=LHmw_Inds+1;

% Get incenters of triangles.
% TR = TriRep(F, V);
% P = TR.incenters;

% vector field
%DirecsVecs=struct('cdata',[],'colormap',[]);
%for j=1:300
%u=us.vf_left{j};
%vATTR=fl.TRs{j};
%figure;
%axis([-1, 1, -1, 1, 0, 1]);
%quiver3(P(:, 1), P(:, 2), P(:, 3), u(:, 1), u(:, 2), u(:, 3), 4, 'k');
%hold on
%trisurf(F, V(:, 1), V(:, 2), V(:, 3), vATTR, 'EdgeColor','none');
%axis equal
%daspect([1, 1, 1]);
%caxis([-45,45]);
%colorbar
%view(115,315); %medial wall view
%view(200,200);
%DirecsVecs(j)=getframe(gcf);
%end
% create videowriter object
%video = VideoWriter('testDirecVecs_rot.avi', 'Uncompressed AVI');
%video.FrameRate = 4;
% open it, plop Direcs in
%open(video)
%writeVideo(video, DirecsVecs);
%close(video);

% color instead of vecs for dir
% compute color space scaling: scaled to last u
%nmax = max(sqrt(sum(u.^2, 2)));
% populate direcs with frames of directionality evolution on sphere
%Direcs=struct('cdata',[],'colormap',[]);
%for j=1:100
% Compute color of projection.
%u=us.vf_left{j};
%c = double(squeeze(computeColour(u(:, 1)/nmax, u(:, 2)/nmax))) ./ 255;
%figure;
%axis([-1, 1, -1, 1, 0, 1]);
%trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
%daspect([1, 1, 1]);
%view(115,315);
%Direcs(j)=getframe(gcf)
%end
% create videowriter object
%video = VideoWriter('testDirec.avi', 'Uncompressed AVI');
%video.FrameRate = 4;
% open it, plop Direcs in
%open(video)
%writeVideo(video, Direcs);
%close(video);

% replicate for original time series

%BOLD=struct('cdata',[],'colormap',[]);
%for j=1:50
%figure;
%axis([-1, 1, -1, 1, 0, 1]);
% verts at this TR
%vATTR=fl.TRs{j};
% apply mask
%vATTR(LHmw_Inds)=-35;
%trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', vATTR, 'EdgeColor', 'none');
%daspect([1, 1, 1]);
%caxis([-45,45]);
%colorbar
% try rotating view
%view(115,315)
%BOLD(j)=getframe(gcf)
%end
% create videowriter object
%video = VideoWriter('testBOLD_rot.avi', 'Uncompressed AVI');
%video.FrameRate = 4;
% open it, plop Direcs in
%open(video)
%writeVideo(video, BOLD);
%close(video);

% right hemi next?

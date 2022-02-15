function us = OpFl_Sph_fs4(subj)

% f should be data in {----}, both in the same order as vx_l (and vx_r) below. 
% adapted from Decomposition of Optical Flow on the Sphere, by Kirisits, Lang, and Scherzer (2014)
% Thank you Kirisits, Lang, and Scherzer!

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
% load in fsaverage5 faces and vertices
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';

%%%%%%%%% HARDCODE FILEPATH TO THESE 
% load in TRs_l and TRs_r
TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_3k.mgh'];
TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_3k.mgh'];
% filepaths to files
TRs_lf=MRIread(TRs_lfp);
TRs_rf=MRIread(TRs_rfp);
% files to data
TRs_l=squeeze(TRs_lf.vol);
TRs_r=squeeze(TRs_rf.vol);
% for surface data
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
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
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
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

% r h 
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
	fr.TRs{TRP}=TRs_r(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
parentfp = '/cbica/projects/hcpd/data/motMasked_contSegs/';
CSIfp = [parentfp subj '/' subj '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'];
CSI = readtable(CSIfp);

% assure that TR count is the same between time series and valid segments txt
SegNum=height(CSI);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI{SegNum,1}+CSI{SegNum,2}-1;
if numTRsVS ~= TR_n
	disp('TRs from Valid Segments txt and cifti do not match. Fix it.')
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow...');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
rmap(roybigbl_cm);

TRPC=1;

% for each continuous segment
for seg=1:SegNum;
	% just to print out count of current segment
	seg
	SegStart=CSI{seg,1};
	SegSpan=CSI{seg,2};
	% get corresponding TRs from aggregate time series
	segTS_l=fl.TRs(SegStart:(SegStart+SegSpan-1));
	segTS_r=fr.TRs(SegStart:(SegStart+SegSpan-1));
	% loop over each TR-Pair: 1 fewer pair than number of TRs
	for TRP=1:(SegSpan-1);
		% print TR pair iter
		TRP
		% Compute decomposition.
		tic;
		% pull out adjacent frames
		u = of(N, faces_l, vx_l, segTS_l{TRP}, segTS_l{TRP+1}, h, alpha, s);
		% throw u into struct
		us.vf_left{TRPC}=u;
		% now right hemi
		u = of(N, faces_l, vx_l, segTS_r{TRP}, segTS_r{TRP+1}, h, alpha, s);
		toc;
		% throw u into struct
		us.vf_right{TRPC}=u;
		% update TR pair counter, which should increase +1 across segments
		TRPC=TRPC+1;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'],'us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualize
% extract medial wall mask
% LHmw=read_label([],'/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/lh.Medial_wall');
% LHmw_Inds=LHmw(:,1);
% might be a +1 in indexing begins at 0 again
% LHmw_Inds=LHmw_Inds+1;

% Get incenters of triangles.
% TR = TriRep(faces_l, vx_l);
% P = TR.incenters;

% vector field
%DirecsVecs=struct('cdata',[],'colormap',[]);
%for j=1:50
%u=us.vf_left{j};
%vATTR=fl.TRs{j};
%figure;
%axis([-1, 1, -1, 1, 0, 1]);
%quiver3(P(:, 1), P(:, 2), P(:, 3), u(:, 1), u(:, 2), u(:, 3), 4, 'k');
%hold on
%trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), vATTR, 'EdgeColor','none');
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


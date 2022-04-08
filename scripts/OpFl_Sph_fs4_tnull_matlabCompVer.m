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
%addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
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

% make a shuffled vector for temporal null
ShufVec=1:TR_n;
% 100 shuffles for now
numShufs=100;

ShufMat=zeros(numShufs,TR_n);
for i=1:numShufs
	ShufMat(i,:)=ShufVec(randperm(TR_n));
end


% initialize TRP counter OUTSIDE OF SHUFFLE LOOP: for plopping u outputs into master struct w/o/r/t their segment
TRPC=1;
% note trp = tr pair

% also initialize output struct outside of shuffle loop
us=struct;

%%%%%%%%%%%%%%%% for each shuffle
for sh=1:numShufs
	tic
	% print shuffle number
	shuffle=sh	
	% reorg time series w/r/t shuffle
	% left hemi
	fl=struct;
	% populate struct
	for TRP=1:TR_n;
		fl.TRs{TRP}=TRs_l(:,ShufMat(sh,TRP));
	end
	
	% r h 
	fr=struct;
	for TRP=1:TR_n;
		fr.TRs{TRP}=TRs_r(:,ShufMat(sh,TRP));
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
	disp('Computing optical flow...');
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
			TRP;
			% Compute decomposition.
			% pull out adjacent frames
			u = of(N, faces_l, vx_l, segTS_l{TRP}, segTS_l{TRP+1}, h, alpha, s);
			% throw u into struct
			us.vf_left{TRPC}=u;
			% now right hemi
			u = of(N, faces_r, vx_r, segTS_r{TRP}, segTS_r{TRP+1}, h, alpha, s);
			% throw u into struct
			us.vf_right{TRPC}=u;
			% update TR pair counter, which should increase +1 across segments
			TRPC=TRPC+1;
		end
	end

toc
%end for each shuffle
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4_t_shuf.mat'],'us','-v7.3')


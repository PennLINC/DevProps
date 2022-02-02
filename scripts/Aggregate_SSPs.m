%%% addpaths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

%%% load in surface data - verts and faces
% for surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

%%%%%%%%%%%%%%%%%%%%%%% Aggregate Subj Loadings
% get subj list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);
% initialize matrix (2 hemis, 17 networks bb)
numFaces=(5120*2)*17;
numSubjs=height(subjs);
FaceTopogMat=zeros(numFaces,numSubjs);
% get indices of where each network should be plopped in to partially vectorized df
Kind_v={};
% starting point for vertex one network 1
Kstart_v=1;
% ending point for network 1
Kend_v=(2562*2);
Kind_v{1}=Kstart_v:Kend_v;
for k=2:17
	% starting point for network k
	Kstart_v=Kend_v+1;
	% ending point for network k
	Kend_v=Kstart_v+(2562*2);
	% K-network indices stored in cell structure
	Kind_v{k}=Kstart_v:Kend_v;
end
% equivalent structure for faces
Kind_f={};
% starting point for vertex one network 1
Kstart_f=1;
% ending point for network 1
Kend_f=(5120*2);
Kind_f{1}=Kstart_f:Kend_f;
for k=2:17
        % starting point for network k
        Kstart_f=Kend_f+1;
        % ending point for network k
        Kend_f=Kstart_f+(5120*2);
        % K-network indices stored in cell structure
        Kind_f{k}=Kstart_f:Kend_f;
end
% for each subj
for s=1:numSubjs
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% for each network
	for k=1:17
	% get indices for this network
	NetInds=Kind_v{k};
	% sep out left and right
	NetIndsL=NetInds(1:2562);
	NetIndsR=NetInds(2563:5124);
	% equivalent for faces
	NetInds_f=Kind_f{k};
        % sep out left and right
        NetIndsL_f=NetInds_f(1:5120);
        NetIndsR_f=NetInds_f(5121:10240);
	% get filepaths
	subjFPL=['/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/resamp_3k_network_L' num2str(k) '.func.gii'];
	subjFPR=['/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/resamp_3k_network_R' num2str(k) '.func.gii'];
	% load in mats
	V_Topog_L=gifti(subjFPL);
	V_Topog_R=gifti(subjFPR);
	V_Topog_L=V_Topog_L.cdata(:,1);
	V_Topog_R=V_Topog_R.cdata(:,1);
	% convert to faces
	F_Topog_L=sum(V_Topog_L(faces_l),2)./3;
	F_Topog_R=sum(V_Topog_R(faces_r),2)./3;
	% plop in facemat!
	FaceTopogMat(NetIndsL_f,s)=F_Topog_L;
	FaceTopogMat(NetIndsR_f,s)=F_Topog_R;
	end
end
%%%%%%%%%%%%%%%%%%%%%

% saveout
save('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjTopog.mat','FaceTopogMat')

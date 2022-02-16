%%% addpaths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

AngD=load([subj '_AngDistMat4.mat']);
file=load([subj '_OpFl_fs4_sigd.mat']);
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
% load in mats
V_L=file.us.sigd_left;
V_R=file.us.sigd_right;
% fs4 faces
facenum=5120;
% tr number
filesize=size(file.us.sigd_left);
% same as below, -1 is tmp
filesize=filesize(2)-1;
%format mats
Facemat_L=zeros(facenum,filesize);
Facemat_R=zeros(facenum,filesize);
% for each TR
for a=1:filesize
	% need +1 because of initial error in code, tmp
	VecVal_L=V_L{a+1};
	VecVal_R=V_R{a+1};
	Facemat_L(:,a)=sum(VecVal_L(faces_l),2)./3;
	Facemat_R(:,a)=sum(VecVal_R(faces_r),2)./3;
end

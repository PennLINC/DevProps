function Gradient_ofPG(subj)
% addpath needed for reading cifti
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% paths needed for reading in surf data
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
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));

%%% read in subject's PG
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];
PG=gifti([childfp sname '_PG_L_10k_rest.func.gii']);

% extract PG values
PG_LH=PG.cdata(:,1);

% get azimuth and elevation from xyz coordinates of sphere
[az,el,r]=cart2sph(vx_l(:,1),vx_l(:,2),vx_l(:,3));
% convert to phi theta: matlab wants us to buy a license for these 4 lines of code?
cos_theta = cos(el).*cos(az);
tan_phi = tan(el) ./ sin(az);
theta = acos(cos_theta);
phi = atan2(tan(el), sin(az));
phi = mod((phi + 2 * pi),(2 * pi));

% to match input style of https://www.mathworks.com/matlabcentral/fileexchange/66482-gradient_sym-v-x-coordinate_system
X=[r,theta,phi];
V=PG_LH;
SphGrad=[diff(V,X(:,1)),1/X(:,1)*diff(V,X(:,2)),1/(X(:,1)*sin(X(:,2)))*diff(V,X(:,3))];




[GxL,GyL]=imgradientxy(masked_PG_gr_L);
[GxR,GyR]=imgradientxy(masked_PG_gr_R);

% save out gradient gradients
dlmwrite([childfp sname '_PG_GxR_BU.csv'],GxR);
dlmwrite([childfp sname '_PG_GyR_BU.csv'],GyR);
dlmwrite([childfp sname '_PG_GxL_BU.csv'],GxL);
dlmwrite([childfp sname '_PG_GyL_BU.csv'],GyL);

% save out low-res gradients for percentile binning
dlmwrite([childfp sname '_PG_lowResFlat_L.csv'],masked_PG_gr_L);
dlmwrite([childfp sname '_PG_lowResFlat_R.csv'],masked_PG_gr_R);

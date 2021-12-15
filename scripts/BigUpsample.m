function BigUpsample(subj,LH,RH)

% addpath needed for reading cifti
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

%%% load in normative x y coordinates for left
FM_l=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.L.flat.32k_fs_LR.surf.gii');
FM_r=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.R.flat.32k_fs_LR.surf.gii');

%%% read in subject's PG for reference
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];
PG=read_cifti([childfp sname '_PG_LR_32k_rest.dscalar.nii']);

% following - https://github.com/coalsont/cifti-matlab/blob/master/cifti_dense_get_surf_map.m
% prepare indices of left hemi
[vertlist1L, ciftilistL, numvertsL] = cifti_dense_get_surf_map(PG.diminfo{1}, 'CORTEX_LEFT');
% prepare indices of right hemi
[vertlist1R, ciftilistR, numvertsR] = cifti_dense_get_surf_map(PG.diminfo{1}, 'CORTEX_RIGHT');

% extract coordinates
xL=double(FM_l.vertices(:,1));
xR=double(FM_r.vertices(:,1));
yL=double(FM_l.vertices(:,2));
yR=double(FM_r.vertices(:,2));

% get indices of where mw verts are
LeftMW=find(sum(FM_l.vertices(:,:)')==0);
RightMW=find(sum(FM_r.vertices(:,:)')==0);

% no negative coordinates
xL=xL+abs(min(xL));
yL=yL+abs(min(yL));
xR=xR+abs(min(xR));
yR=yR+abs(min(yR));

% going to need to expand these to make seperable bins IN ACCORDANCE WITH PROC. POWER AVAILABLE
SxL=xL*.2; % example: xL*10 for 10x resolution
SxR=xR*.2;
SyL=yL*.2;
SyR=yR*.2;

%%% just need to sub in x,y coords for 64*89 as first and second arg to griddata, grid data in same order, and xL yL and 3rd and 4th

% coords of data
[Lrow,Lcol]=find(~isnan(LH));
[Rrow,Rcol]=find(~isnan(RH));

% only valid pixels
%for P=1:length(Lrow)
%LH_v(P)=LH(Lrow(P),Lcol(P));
%end

% mirror for right


% rescale image to size of cifti flatmap
LH_big = imresize(LH,[max(yL),max(xL)]);

% coords of resamp data
[Lrow,Lcol]=find(~isnan(LH_big));

% same loop to vectoriz big map
LH_bigV=[];
for P=1:length(Lrow)
LH_bigV(P)=LH_big(Lrow(P),Lcol(P));
end

LH_flatmap=griddata(double(Lrow),double(Lcol),double(LH_bigV),double(yL),double(xL));

% vert ordering
LH_ordered=1:max(vertlist1L);
LH_ordered(vertlist1L)=LH_flatmap(vertlist1L);
LH_ordered(LeftMW)=0;

test=[];
% * 26 just to humor cifti
for D=1:26
test(:,D)=LH_ordered;
end
a=cifti_dense_replace_surfdata(PG,test,'CORTEX_LEFT',1)
write_cifti(a,'testcif.dscalar.nii')  

%%%% STOP IT!
%%% get flatmap and masks in order
% need to use combined index where xL yL xR yR are not 0 and cifti index is not 0, first step is filtering by vertlist from ciftiinfo
% this extracts the shape of the flatmap in the form of x and y coordinates - vertices of all coordinates comprising the flatmaps
xLPartialFilt=xL(vertlist1L);
yLPartialFilt=yL(vertlist1L);
xRPartialFilt=xR(vertlist1R);
yRPartialFilt=yR(vertlist1R);
% no negative coordinates
xLPartialFilt=xLPartialFilt+abs(min(xLPartialFilt));
yLPartialFilt=yLPartialFilt+abs(min(yLPartialFilt));
xrRPartialFilt=xRPartialFilt+abs(min(xRPartialFilt));
yRPartialFilt=yRPartialFilt+abs(min(yRPartialFilt));
% size of "grid"
s_gridX_L=max(xL)-min(xL);
s_gridX_R=max(xR)-min(xR);
s_gridY_L=max(yL)-min(yL);
s_gridY_R=max(yR)-min(yR);
% Xq and Yq will be equally spaced integer grid values
XqL=1:s_gridX_L;
YqL=1:s_gridY_L;
XqR=1:s_gridX_R;
YqR=1:s_gridY_R;
% convert to grid
%[xL,yL]=meshgrid(XqL,YqL);
%[xR,yR]=meshgrid(XqR,YqR);
% extract PG values
PG_LH=PG.cdata(ciftilistL,1);
% extract right
PG_RH=PG.cdata(ciftilistR,1);


% C is 64*89 in this toy example
% source of transpose is matlab and python naturally don't represent matrices the same way: one is transposed from the other 
% JUST FILL IN 0S for ALL MW VERTS IN 32K BUT NOT 29K - USE 
PG_gr_L = griddata(double(xL),double(yL),double(LH),double(xLPartialFilt),double(yLPartialFilt));

cifti_dense_replace_surfdata(PG,PG_gr_L,'CORTEX_LEFT',1)

PG_gr_L = griddata(double(xL),double(yL),double(C),double(xLPartialFilt),double(yLPartialFilt));
PG_gr_R = griddata(double(xRPartialFilt),double(yRPartialFilt),double(PG_RH),double(xR),double(yR));


% would try https://github.com/coalsont/cifti-matlab/blob/master/cifti_dense_replace_surfdata.m first

% lifted from abcdfnets pipeline
addpath(genpath('/cbica/projects/abcdfnets/scripts/cifti-matlab/'));
% cifti to replace cdata in
HP=cifti_read('/cbica/projects/abcdfnets/results/SingleParcellation/RobustInitialization_Cifti_Surf/robust_initVHP.dscalar.nii');
HP.cdata(1:59412)=sbj_AtlasLabel_NoMedialWall;
outputfile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/' subj '_Parcel.dscalar.nii'];
cifti_write(HP,outputfile)



% save out gradient gradients
dlmwrite([childfp sname '_PG_GxR_BU.csv'],GxR);
dlmwrite([childfp sname '_PG_GyR_BU.csv'],GyR);
dlmwrite([childfp sname '_PG_GxL_BU.csv'],GxL);
dlmwrite([childfp sname '_PG_GyL_BU.csv'],GyL);

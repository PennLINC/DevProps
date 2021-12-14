function permut_PGs(subj)
% create the permutation (null) maps of hierarchy directionality

% spin the gradient 1000 times
spin_pg(subj)

% set FP
sname=char(subj);
parentfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];

%%% load in resources needed for resampling
% addpath needed for reading cifti
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% and path needed for opflow
addpath(genpath('/cbica/projects/abcdfnets/scripts/NeuroPattToolbox'));
% load in normative x y coordinates for left
FM_l=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.L.flat.32k_fs_LR.surf.gii');
FM_r=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.R.flat.32k_fs_LR.surf.gii');
% extract coordinates
xL=double(FM_l.vertices(:,1));
xR=double(FM_r.vertices(:,1));
yL=double(FM_l.vertices(:,2));
yR=double(FM_r.vertices(:,2));
% make downsampled grid for interpolation
xL=xL*.2; % example: xL*10 for 10x resolution
xR=xR*.2;
yL=yL*.2;
yR=yR*.2;
% original PG for reference
PG=read_cifti([parentfp sname '_PG_LR_32k_rest.dscalar.nii']);
% prepare indices of left hemi
[vertlist1L, ciftilistL, numvertsL] = cifti_dense_get_surf_map(PG.diminfo{1}, 'CORTEX_LEFT');
% prepare indices of right hemi
[vertlist1R, ciftilistR, numvertsR] = cifti_dense_get_surf_map(PG.diminfo{1}, 'CORTEX_RIGHT');
xLPartialFilt=xL(vertlist1L);
yLPartialFilt=yL(vertlist1L);
xRPartialFilt=xR(vertlist1R);
yRPartialFilt=yR(vertlist1R);
% no negative coordinates
xLPartialFilt=xLPartialFilt+abs(min(xLPartialFilt));
yLPartialFilt=yLPartialFilt+abs(min(yLPartialFilt));
xRPartialFilt=xRPartialFilt+abs(min(xRPartialFilt));
yRPartialFilt=yRPartialFilt+abs(min(yRPartialFilt));
% size of "grid"
s_gridX_L=max(xL)-min(xL);
s_gridX_R=max(xR)-min(xR);
s_gridY_L=max(yL)-min(yL);
s_gridY_R=max(yR)-min(yR);
% get border of flatmap in a sized rectangle
vqBound_L = boundary(double(xLPartialFilt),double(yLPartialFilt)); 
vqBound_R = boundary(double(xRPartialFilt),double(yRPartialFilt));	
% Xq and Yq will be equally spaced integer grid values
XqL=1:s_gridX_L;
YqL=1:s_gridY_L;
XqR=1:s_gridX_R;
YqR=1:s_gridY_R;
% convert to grid
[xL,yL]=meshgrid(XqL,YqL);
[xR,yR]=meshgrid(XqR,YqR);
% create mask for vertices within boundaries of shape (within rectangle)
bwL = poly2mask(double(xLPartialFilt(vqBound_L)),double(yLPartialFilt(vqBound_L)),double(max(max(yL))),double(max(max(xL))));
bwR = poly2mask(double(xLPartialFilt(vqBound_R)),double(yLPartialFilt(vqBound_R)),double(max(max(yR))),double(max(max(xR))));	
% convert to NaN's instead of 0
bwNL=double(bwL);
bwNL(bwNL==0)=NaN;
bwNR=double(bwR);
bwNR(bwNR==0)=NaN;

% initialize struct: instead of a 5d array, setting Left and Right hemis to have sep. X and Y gradient fields within struct
SpunPG_Gs=struct();
SpunPG_Gs.Lefts={};
SpunPG_Gs.Rights={};
% note left and right flatmap have slightly different x-axis length
SpunPG_Gs.Lefts.Xs=zeros(1000,64,89);
SpunPG_Gs.Rights.Xs=zeros(1000,68,89);
SpunPG_Gs.Lefts.Ys=zeros(1000,64,89);
SpunPG_Gs.Rights.Ys=zeros(1000,68,89);

% For each spin, project to high-res, get gradient gradient, insert into struct, delete high-res map. 
for s=1:1000
	% just for some prinout to keep tabs
	iteration=s

	% convert data vector to brainmap (fsaverage5 space)
	applySpinsCmd=['/cbica/projects/pinesParcels/miniconda3/envs/mv_preds/bin/python applyPGSpins.py ' sname ' ' num2str(s) ];
	system(applySpinsCmd);
	%pause(xxx)

	% project to high-res
	UScmd=['./upsample_PG_spun.sh ' sname];
	system(UScmd);
	%pause(xxx)

	% load in upsample
	usFN=[parentfp sname '_PG_LR_spun_32k.dscalar.nii'];
	spunion=read_cifti(usFN);

	% extract PG values
	PG_LH=spunion.cdata(ciftilistL,1);
	% extract right
	PG_RH=spunion.cdata(ciftilistR,1);

	% downsamplesample to flatmap
	PG_gr_L = griddata(double(xLPartialFilt),double(yLPartialFilt),double(PG_LH),double(xL),double(yL));
	PG_gr_R = griddata(double(xRPartialFilt),double(yRPartialFilt),double(PG_RH),double(xR),double(yR));
	
	% get the gradient of the PG (I know, language sucks)
	[GxL,GyL]=imgradientxy(PG_gr_L);
	[GxR,GyR]=imgradientxy(PG_gr_R);

	% insert resample into struct
	SpunPG_Gs.Left.Xs(s,:,:)=GxL;
	SpunPG_Gs.Left.Ys(s,:,:)=GyL;
        SpunPG_Gs.Right.Xs(s,:,:)=GxR;
        SpunPG_Gs.Right.Ys(s,:,:)=GyR;
end

outFN=[parentfp sname '_spunFlats.mat'];
save(outFN,'SpunPG_Gs')


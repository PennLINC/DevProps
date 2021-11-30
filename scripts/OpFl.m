function OpFl(subj)
% run optical flow analyses on this subject, on each segment of continuous TRs sep. Only on resting state for now. 
% main to-do's:
% find way to mask empty cells
% hopefully find a way to speed it up
% find a way to combine SVD results from each run 
% re=semicolon unsemicolon'ed lines


% addpath needed for reading ciftis and opflow
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

%%% load in normative x y coordinates for left
FM_l=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.L.flat.32k_fs_LR.surf.gii');
FM_r=gifti('/cbica/projects/pinesParcels/data/Surfs/S900.R.flat.32k_fs_LR.surf.gii');
% extract coordinates
xL=double(FM_l.vertices(:,1));
xR=double(FM_r.vertices(:,1));
yL=double(FM_l.vertices(:,2));
yR=double(FM_r.vertices(:,2));
% going to need to expand these to make seperable bins IN ACCORDANCE WITH PROC. POWER AVAILABLE
xL=xL*.1; % example: xL*10 for 10x resolution
xR=xR*.1;
yL=yL*.1;
yR=yR*.1;

%%% read in subject's clean TS, starting with rest
sname=char(subj);
parentfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' sname '/'];
ts=read_cifti([parentfp sname '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii']);

% following - https://github.com/coalsont/cifti-matlab/blob/master/cifti_dense_get_surf_map.m
% prepare indices of left hemi
[vertlist1L, ciftilistL, numvertsL] = cifti_dense_get_surf_map(ts.diminfo{1}, 'CORTEX_LEFT');

% prepare indices of right hemi
[vertlist1R, ciftilistR, numvertsR] = cifti_dense_get_surf_map(ts.diminfo{1}, 'CORTEX_RIGHT');

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
%bwNL(bwNL==0)=NaN;
bwNR=double(bwR);
%bwNR(bwNR==0)=NaN;

%%% get time series in order
% extract left hemi
ts_LH=zeros(numvertsL,1,'single');
ts_LH=ts.cdata(ciftilistL,:);
% extract right
ts_RH=zeros(numvertsR,1,'single');
ts_RH=ts.cdata(ciftilistR,:);
% read in continuous segment denotation
segmentFN=[parentfp sname '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'];
segments=dlmread(segmentFN);
num_segs=length(segments);

%%% set opflow params
params = setNeuroPattParams(1.25); 
params = setNeuroPattParams(params,'zscoreChannels', 1, 1.25);
% params = setNeuroPattParams(params,'subtractBaseline', 0, 1.25); 
% filter does not seem to be stopping 
% params = setNeuroPattParams(params,'filterData',logical(0), 1.25);
% if you can't beat 'em, join em
params = setNeuroPattParams(params,'morletCfreq', .05, 1.25);
params = setNeuroPattParams(params,'opBeta', 10, 1.25);
params = setNeuroPattParams(params,'planeWaveThreshold', 0.7, 1.25);
params = setNeuroPattParams(params,'synchronyThreshold', 0.7, 1.25);
params = setNeuroPattParams(params,'minDurationSecs', 10, 1.25);
params = setNeuroPattParams(params,'maxTimeGapSecs', 5, 1.25);
params = setNeuroPattParams(params,'maxDisplacement', 1, 1.25);
params = setNeuroPattParams(params,'minCritRadius', 1, 1.25);
% setting opAlpha to 1 to account for increased spatial domain/global patterning
params.opAlpha =1.5;
% downsampling the temporal domain by 5x hurts though, they say default is 1
params = setNeuroPattParams(params,'downsampleScale', 1, 1.25);
% initialize megastruct for results from each segment
MegaStruct=struct();
MegaStruct.P_Left={};
MegaStruct.P_Right={};
Mega.TRsInSeg={};
MegaStruct.Vf_Left={};
MegaStruct.Vf_Right={};
%%% for each contin. segment, run opflow
for S=1:num_segs
	% get length of segment in TR
	ts_seg_length=segments(S,2);
	% index into master time series to grab this segment
	ts_segL=ts_LH(:,(segments(S,1):((segments(S,1)+segments(S,2))-1)));
	ts_segR=ts_RH(:,(segments(S,1):((segments(S,1)+segments(S,2))-1)));
	% initialize 3d (x,y,time) flatmaps for both hemis
	fMap_ts_L_grid=zeros(length(YqL),length(XqL),ts_seg_length);
	fMap_ts_R_grid=zeros(length(YqR),length(XqR),ts_seg_length);
	% insert "gridded" TRs into matrix
	for T=1:ts_seg_length
		% extract this TR
		TRVals_L=ts_segL(:,T);
		TRVals_R=ts_segR(:,T);
		% Interp. onto grid CHECK MESHGRID TO LOOK FOR SOURCE OF TRANSPOSE HERE
		vqL = griddata(double(xLPartialFilt),double(yLPartialFilt),double(TRVals_L),double(xL),double(yL));
		vqR = griddata(double(xRPartialFilt),double(yRPartialFilt),double(TRVals_R),double(xR),double(yR));
		% use mask to select flatmap patch
		masked_vqL=(vqL).*(bwNL);
		masked_vqR=(vqR).*(bwNR);
		% insert into broader DF
		fMap_ts_L_grid(:,:,T)=masked_vqL;
		fMap_ts_R_grid(:,:,T)=masked_vqR;
	end	
	% small patch test
	% resultsL = mainProcessingWithOutput(fMap_ts_L_grid(100:140,100:140,:), 1.25, params);
	% run opflow
	resultsL = mainProcessingWithOutput(fMap_ts_L_grid, 1.25, params);
	resultsR = mainProcessingWithOutput(fMap_ts_R_grid, 1.25, params);
	% extract opflow results of interest
	MegaStruct.P_Left{S}=resultsL.patterns;
	MegaStruct.P_Right{S}=resultsR.patterns;
        MegaStruct.Vf_Left{S}=resultsL.velocityFields;
        MegaStruct.Vf_Right{S}=resultsR.velocityFields;
	MegaStruct.TRsInSeg{S}=resultsL.nTimeSteps+1
 end

% make a subject-level directory
s_levelDir=['/cbica/projects/pinesParcels/results/OpFl_output/' sname];
mkdircommand=['mkdir ' s_levelDir];
system(mkdircommand)

% save aggregated output into waveoutput
fn=['/cbica/projects/pinesParcels/results/OpFl_output/' sname '/OpFlowResults.mat'];
save(fn,'MegaStruct')


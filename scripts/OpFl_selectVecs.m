function OpFl_selectVecs(subj)

% pull out a few choice points for circular distributions of vectors

% addpath needed for reading cifti
addpath(genpath('/cbica/projects/abcdfnets/scripts/code_nmf_cifti/tool_folder'));
% and path needed for opflow
addpath(genpath('/cbica/projects/abcdfnets/scripts/NeuroPattToolbox'));

%%% load in normative x y coordinates for left
FM_l=gifti('/cbica/projects/abcdfnets/scripts/normative_surfs/S900.L.flat.32k_fs_LR.surf.gii');
FM_r=gifti('/cbica/projects/abcdfnets/scripts/normative_surfs/S900.R.flat.32k_fs_LR.surf.gii');
% extract coordinates
xL=double(FM_l.vertices(:,1));
xR=double(FM_r.vertices(:,1));
yL=double(FM_l.vertices(:,2));
yR=double(FM_r.vertices(:,2));
% going to need to expand these to make seperable bins IN ACCORDANCE WITH PROC. POWER AVAILABLE
xL=xL*.2; % example: xL*10 for 10x resolution
xR=xR*.2;
yL=yL*.2;
yR=yR*.2;

%%% read in subject's clean TS, starting with rest
sname=char(subj);
parentfp=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];

% load aggregated output from waveoutput
fn=['/cbica/projects/abcdfnets/results/wave_output/' sname '/' 'OpFlowResults.mat'];
results=load(fn);

% get number of runs
numRuns=length(results.MegaStruct.TRsInSeg)

% point 1, mid range on PG
% point1=[30 30];
% point 2 is high on PG
% point2=[60 71];
% point 3 is low on PG
% point3=[46 56];

% initialize vector of vector angles 
p1Vec=[];
p2Vec=[];
p3Vec=[];

% for each run, extract vectors assigned to each point and aggregate
for R=1:numRuns
	Vf_L=results.MegaStruct.Vf_Left{R};
	p1Vec=[p1Vec ; squeeze(Vf_L(30,30,:))];
	p2Vec=[p2Vec ; squeeze(Vf_L(60,71,:))];
	p3Vec=[p3Vec ; squeeze(Vf_L(46,56,:))];
end

% save out vectors of angles
p1FN=['/cbica/projects/abcdfnets/results/wave_output/' sname '/' sname '_p1Vec.csv'];
dlmwrite(p1FN,p1Vec);
p2FN=['/cbica/projects/abcdfnets/results/wave_output/' sname '/' sname '_p2Vec.csv'];
dlmwrite(p2FN,p2Vec);
p3FN=['/cbica/projects/abcdfnets/results/wave_output/' sname '/' sname '_p3Vec.csv'];
dlmwrite(p3FN,p3Vec);


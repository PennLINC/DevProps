% add paths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% read in subj list
Subjs=readtable('~/PWs/rs_subs.csv');
% set parent FP for shorter FPs in loop
parentfp = '/cbica/projects/pinesParcels/results/PWs/Proced/';
% initialize group average matrices
groL=zeros(4589);
groR=zeros(4595);
% for each subj
for s=1:height(Subjs)
	subj=table2array(Subjs{s,2});
	% load in CFC mat
	CFC_fp=[parentfp subj '/' subj '_CircFC.mat'];
	CFC=load(CFC_fp);
	CFC_L=CFC.adjMats.L;
	CFC_R=CFC.adjMats.R;
	groL=groL+CFC_L;
	groR=groR+CFC_R;
	s
	CFC_L(3230,1)
end
% average
groL=groL./height(Subjs);
groR=groR./height(Subjs);

% save out FC mat as .mat
outFN=['/cbica/projects/pinesParcels/results/PWs/Gro_CFC_L.mat'];
save(outFN,'groL')
outFN=['/cbica/projects/pinesParcels/results/PWs/Gro_CFC_R.mat'];
save(outFN,'groR')


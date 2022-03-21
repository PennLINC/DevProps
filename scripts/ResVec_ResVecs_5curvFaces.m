% extract column 8, mean resultant vector, from each angular distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for circle stat toolbox
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% load in subjects list
subjs=readtable('~/PWs/rs_subs.csv')

% initialize output dfs
FaceMatL=zeros(height(subjs),20480);
FaceMatR=zeros(height(subjs),20480);
FaceVecL=zeros(1,20480);
FaceVecR=zeros(1,20480);

% for each subj
numSubjs=height(subjs);
for s=1:numSubjs
        % get subject ID
	subjcell=table2array(subjs(s,2));
        subj=subjcell{:}
	% get resultant vecs for this subj
	FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_curv5.mat'];
	FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_curv5.mat'];
	FileL=load(FP_L);
	FileR=load(FP_R);
	% extract resvecs - 8th column
	BUTDL=FileL.OutDf_L;
	BUTDR=FileR.OutDf_R;
	% always a process with matlab
	FaceMatL(s,:)=cell2mat(BUTDL(:,8));
	FaceMatR(s,:)=cell2mat(BUTDR(:,8));
end

% get mean resultant vector over subjs
% for each face
for f=1:20480;
	FaceVecL(f)=circ_mean(FaceMatL(:,f));
end

% right hemi
for f=1:20480;
	FaceVecR(f)=circ_mean(FaceMatR(:,f));
end

% print into facewise
Vis_FaceVec(FaceVecL,FaceVecR,'GrandMeanMapCurv.png')

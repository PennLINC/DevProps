% extract column 8, mean resultant vector, from each angular distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for circle stat toolbox
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% load in subjects list
subjs=readtable('~/PWs/rs_subs.csv')

% initialize output dfs
FaceMatL=zeros(height(subjs),5120);
FaceMatR=zeros(height(subjs),5120);
FaceVecL=zeros(1,5120);
FaceVecR=zeros(1,5120);

% for each subj
numSubjs=height(subjs);
for s=1:numSubjs
        % get subject ID
	subjcell=table2array(subjs(s,2));
        subj=subjcell{:}
	% get resultant vecs for this subj
	FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_curv4.mat'];
	FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_curv4.mat'];
	% if exists
	if isfile([FP_L]);
		FileL=load(FP_L);
		FileR=load(FP_R);
		% extract resvecs - 8th column
		BUTDL=FileL.OutDf_L;
		BUTDR=FileR.OutDf_R;
		% replace empties with 0
		emptyIndex = cellfun('isempty', BUTDL(:,8));     % Find indices of empty cells
		BUTDL(emptyIndex,8) = {0};
		emptyIndex = cellfun('isempty', BUTDR(:,8));     % Find indices of empty cells
	        BUTDR(emptyIndex,8) = {0};
		% always a process with matlab
		FaceMatL(s,:)=cell2mat(BUTDL(:,8));
		FaceMatR(s,:)=cell2mat(BUTDR(:,8));
	end
end

% get mean resultant vector over subjs
% for each face
for f=1:5120;
	FaceVecL(f)=circ_mean(FaceMatL(:,f));
end

% right hemi
for f=1:5120;
	FaceVecR(f)=circ_mean(FaceMatR(:,f));
end

% print into facewise
Vis_FaceVec(FaceVecL,FaceVecR,'GrandMeanMap4Curv.png')

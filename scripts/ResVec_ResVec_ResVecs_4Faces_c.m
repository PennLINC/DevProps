% resultant vectors to collapse across TRs, subjects, AND faces. Prinout in python friendly format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for circle stat toolbox
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% load in subjects list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G300_cTRs.txt','ReadVariableNames',false);
% initialize output dfs
MasterMat=zeros(height(subjs),(4851+4842));
MasterVector=zeros(1,(height(subjs)*(4851+4842)));

% for each subj
numSubjs=height(subjs);
for s=1:numSubjs
        % get subject ID
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% get resultant vecs for this subj
	FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_c.mat'];
	FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_c.mat'];
	FileL=load(FP_L);
	FileR=load(FP_R);
	% extract resvecs - 8th column
	BUTDL=FileL.OutDf_L;
	BUTDR=FileR.OutDf_R;
	% always a process with matlab
	ResVecsL=cell2mat(BUTDL(:,8));
	ResVecsR=cell2mat(BUTDR(:,8));
	% plop em in DF
	ResVecs=vertcat(ResVecsL,ResVecsR);
	MasterMat(s,:)=ResVecs;
end

% convert mastermat to vector
MMsize = size(MasterMat);
MasterMatElementsNumber = MMsize(1)*MMsize(2);
MasterVector=reshape(MasterMat,1,MasterMatElementsNumber);

% print grand mean mean
GMM_ResVec=circ_mean(MasterVector)

% saveout in python-friendly form.
writetable(table(MasterVector),'~/results/PWs/masterResVecVec_c.csv')

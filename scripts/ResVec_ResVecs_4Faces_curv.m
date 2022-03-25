% extract column 8, mean resultant vector, from each angular distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for circle stat toolbox
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% load in subjects list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);

% initialize output dfs
FaceMatL=zeros(height(subjs),5120);
FaceMatR=zeros(height(subjs),5120);
FaceVecL=zeros(1,5120);
FaceVecR=zeros(1,5120);
stdFaceVecL=zeros(1,5120);
stdFaceVecR=zeros(1,5120);

% for each subj
numSubjs=height(subjs);
for s=1:numSubjs
        % get subject ID
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% get resultant vecs for this subj
	FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_curv4.mat'];
	FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_curv4.mat'];
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
for f=1:5120;
	FaceVecL(f)=circ_mean(FaceMatL(:,f));
	stdFaceVecL(f)=circ_std(FaceMatL(:,f));	
end

% right hemi
for f=1:5120;
	FaceVecR(f)=circ_mean(FaceMatR(:,f));
	stdFaceVecR(f)=circ_std(FaceMatR(:,f));	
end

% print into facewise
Vis_FaceVec_circ(FaceVecL,FaceVecR,'GrandMeanMap.png')
Vis_FaceVec_lin(stdFaceVecL,stdFaceVecR,'GrandsdMap.png')
% combine mean and SD's into a table

% save out distributions
writetable(table(FaceVecL,FaceVecR,stdFaceVecL,stdFaceVecR),'~/results/MeanSDs_MeanRelativeAngle_PGG_L.csv')



%%% Single subject portion

% manually set subj: can't upload to github with subj id
subj='';

% same code as above
FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_curv4.mat'];
FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_curv4.mat'];
FileL=load(FP_L);
FileR=load(FP_R);
% extract resvecs - 8th column
BUTDL=FileL.OutDf_L;
BUTDR=FileR.OutDf_R;
% always a process with matlab
FaceVecL=cell2mat(BUTDL(:,8));
FaceVecR=cell2mat(BUTDR(:,8));
stdFaceVecL=cell2mat(BUTDL(:,9));
stdFaceVecR=cell2mat(BUTDR(:,9));
Vis_FaceVec_circ(FaceVecL,FaceVecR,'SubjMeanMap_curv.png')
Vis_FaceVec_lin(stdFaceVecL,stdFaceVecR,'SubjsdMap_curv.png')

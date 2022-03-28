% extract column 8, mean resultant vector, from each angular distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for circle stat toolbox
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% load in subjects list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);

% initialize output dfs
FaceMatL=zeros(height(subjs),4589);
FaceMatR=zeros(height(subjs),4595);
FaceVecL=zeros(1,4589);
FaceVecR=zeros(1,4595);
stdFaceVecL=zeros(1,4589);
stdFaceVecR=zeros(1,4595);
% Arbitrary Direction Mean Vec too
ArbDFaceMatL=zeros(height(subjs),4589);
ArbDFaceMatR=zeros(height(subjs),4595);
ArbDFaceVecL=zeros(1,4589);
ArbDFaceVecR=zeros(1,4595);
% for each subj
numSubjs=height(subjs);
for s=1:numSubjs
        % get subject ID
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% get resultant vecs for this subj
	FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
	FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
	FileL=load(FP_L);
	FileR=load(FP_R);
	% extract resvecs - 8th column
	BUTDL=FileL.OutDf_L;
	BUTDR=FileR.OutDf_R;
	% always a process with matlab
	FaceMatL(s,:)=cell2mat(BUTDL(:,8));
	FaceMatR(s,:)=cell2mat(BUTDR(:,8));
	% mean direction relative to spherical north (y/elevation = 1, x/azimuth=0)
	ArbDFaceMatL(s,:)=cell2mat(BUTDL(:,10));
	ArbDFaceMatR(s,:)=cell2mat(BUTDR(:,10));
end

% get mean resultant vector over subjs
% for each face
for f=1:4589;
	FaceVecL(f)=circ_mean(FaceMatL(:,f));
	stdFaceVecL(f)=circ_std(FaceMatL(:,f));	
	ArbDFaceVecL(f)=circ_mean(ArbDFaceMatL(:,f));
end

% right hemi
for f=1:4595;
	FaceVecR(f)=circ_mean(FaceMatR(:,f));
	stdFaceVecR(f)=circ_std(FaceMatR(:,f));	
	ArbDFaceVecR(f)=circ_mean(ArbDFaceMatR(:,f));
end

% print into facewise
Vis_FaceVec_circ(FaceVecL,FaceVecR,'GrandMeanMap.png')
Vis_FaceVec_lin(stdFaceVecL,stdFaceVecR,'GrandsdMap.png')
Vis_FaceVec_circ(ArbDFaceVecL,ArbDFaceVecR,'ArbDMeanMap.png')

% save out distributions
writetable(table(FaceVecL,FaceVecR,stdFaceVecL,stdFaceVecR,ArbDFaceVecL,ArbDFaceVecR),'~/results/MeanSDs_MeanRelativeAngle_PGG.csv')


%%% Single subject portion

% manually set subj: can't upload to github with subj id
%subj='';

% same code as above
FP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
FP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
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
Vis_FaceVec_circ(wrapToPi(FaceVecL),wrapToPi(FaceVecR),'SubjMeanMap.png')
Vis_FaceVec_lin(stdFaceVecL,stdFaceVecR,'SubjsdMap.png')
% arbitrary dir circ faceplot
Vis_FaceVec_circ(cell2mat(BUTDL(:,10)),cell2mat(BUTDR(:,10)),'SubjMeanMap_arbDInfl.png')

%%% addpaths

%%%%%%%%%%%%%%%%%%%%%%% Aggregate Subj Loadings
% get subj list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);
% initialize matrix (2 hemis, 4 components per face)
numFaces=(5120*2)*5;
numSubjs=height(subjs);
FaceVecCompMat=zeros(numFaces,numSubjs);
% create indices dictating where to plop each vector component in the long matrix

BULH=1:5120;
BULV=5121:10240;
BURH=10241:15360;
BURV=15361:20480;
TDLH=20481:25600;
TDLV=25601:30720;
TDRH=30721:35840;
TDRV=35841:40960;
PropBUL=40961:46080;
PropBUR=46081:51200;

% for each subj
for s=1:numSubjs
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% get filepaths
	subjFPL=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
	subjFPR=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
	% load in mats
	BUTDLdf=load(subjFPL);
	BUTDRdf=load(subjFPR);
	% extract structs
	BUTDL=BUTDLdf.OutDf_L;
	BUTDR=BUTDRdf.OutDf_R;
	% replace empties with 0s
	emptiesL=cellfun('isempty',BUTDL);
	BUTDL(emptiesL)={0};
	emptiesR=cellfun('isempty',BUTDR);
        BUTDR(emptiesR)={0};
	% grab bottom up components
	FaceVecCompMat(BULH,s)=cell2mat(BUTDL(:,3));
	FaceVecCompMat(BURH,s)=cell2mat(BUTDR(:,3));
	FaceVecCompMat(BULV,s)=cell2mat(BUTDL(:,2));
        FaceVecCompMat(BURV,s)=cell2mat(BUTDR(:,2));
	% and top down
	FaceVecCompMat(TDLH,s)=cell2mat(BUTDL(:,5));
        FaceVecCompMat(TDRH,s)=cell2mat(BUTDR(:,5));
        FaceVecCompMat(TDLV,s)=cell2mat(BUTDL(:,4));
        FaceVecCompMat(TDRV,s)=cell2mat(BUTDR(:,4));
	% and prop
	FaceVecCompMat(PropBUL,s)=cell2mat(BUTDL(:,1));
        FaceVecCompMat(PropBUR,s)=cell2mat(BUTDR(:,1));
end
%%%%%%%%%%%%%%%%%%%%%

% saveout
save('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjVecComps.mat','FaceVecCompMat')

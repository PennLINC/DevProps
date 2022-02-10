% use PLS weights and PCA decomp. to reconstruct PLS weights in native feature space
% load pls weights
SSPWeightsM=load('/cbica/projects/pinesParcels/results/PWs/SSP_Weights.mat');
VecWeightsM=load('/cbica/projects/pinesParcels/results/PWs/Vec_Weights.mat');
SSPWeights=SSPWeightsM.SSP_Weights;
VecWeights=VecWeightsM.Vec_Weights;
% load PCA decomps
SSP_PCAM=load('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjTopogPCA.mat');
SSP_PCA=SSP_PCAM.SSP_PC_struct;
Vec_PCAM=load('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjVecsPCA.mat');
Vec_PCA=Vec_PCAM.Vecs_PC_struct;

% adapt to number of final selected PCs
numPCs=length(SSPWeights)

% initialize feature weight matrix for SSPs and vectors
SSP_featWeights=zeros(numPCs,174080);
Vec_featWeights=zeros(numPCs,40960);

% reconstruct into original feature space, by adding the absolute feature weight of each PC in accordance with its PC loading

% PLS component number
PLcN=1

% for each principal component
for P=1:numPCs
	% extract PC weights from this PLS component
	SSP_PCws=SSPWeights(P,PLcN);
	% average over all iterations
	%av_SSP_PCws=mean(mean((SSP_PCws)));
	% same for vectors
	Vec_PCws=VecWeights(P,PLcN);
	% average over all iterations
	%av_Vec_PCws=mean(mean((Vec_PCws)));
	% multiply by PC weights to get into feature space
	% pc weights in abs.
	SSP_featWeights(P,:)=(SSP_PCws*SSP_PCA.coeff(:,P));
	Vec_featWeights(P,:)=(Vec_PCws*Vec_PCA.coeff(:,P));
	P
end

% get mean weight for each native feature across PCs
mSSP_featWeights=mean(SSP_featWeights);
mVec_featWeights=mean(Vec_featWeights);

%%%% index out for printing % - copied from PCA scripts
% Vecs
BULH=1:5120;
BULV=5121:10240;
BURH=10241:15360;
BURV=15361:20480;
TDLH=20481:25600;
TDLV=25601:30720;
TDRH=30721:35840;
TDRV=35841:40960;

% Nets - build the indices for extraction later
% same indices as before
% equivalent structure for faces
Kind_f={};
% starting point for vertex one network 1
Kstart_f=1;
% ending point for network 1
Kend_f=(5120*2);
Kind_f{1}=Kstart_f:Kend_f;
for k=2:17
        % starting point for network k
        Kstart_f=Kend_f+1;
        % ending point for network k
        Kend_f=Kstart_f+(5120*2)-1;
        % K-network indices stored in cell structure
        Kind_f{k}=Kstart_f:Kend_f;
end

% extract, combine, TD BU horz and vert
BUL=zeros(1,5120);
BUR=zeros(1,5120);
TDL=zeros(1,5120);
TDR=zeros(1,5120);

% just some annoying acronym overload here, no way around it
BUL(:)=mVec_featWeights(BULH)+mVec_featWeights(BULV);
BUR(:)=mVec_featWeights(BURH)+mVec_featWeights(BURV);
TDL(:)=mVec_featWeights(TDLH)+mVec_featWeights(TDLV);
TDR(:)=mVec_featWeights(TDRH)+mVec_featWeights(TDRV);

Vis_FaceVec(BUL,BUR,'BU_PLSc1_ws.png')
Vis_FaceVec(TDL,TDR,'TD_PLSc1_ws.png')

% merged weight vector
MergedL=zeros(1,5120);
MergedR=zeros(1,5120);

% extract out each net and vis 
for k=1:17
        % get indices belonging to this network
        Inds=Kind_f{k};
        % extract the left face for this network
        NetIndsL=Inds(1:5120);
        FaceVecL=mSSP_featWeights(NetIndsL);
        % extract the right face for this network
        NetIndsR=Inds(5121:10240);
        FaceVecR=mSSP_featWeights(NetIndsR);
        % print mean coefficient loadings onto this network
        mean(abs(FaceVecR))
        % create png on surface
        Fn=['~/results/PWs/PLSc1_Net' num2str(k) '_posNeg.png'];
        Vis_FaceVec(FaceVecL,FaceVecR,Fn)
	MergedL=MergedL+FaceVecL;
	MergedR=MergedR+FaceVecR;
end

Vis_FaceVec(MergedL,MergedR,'~/results/PWs/PLSc1_mergedNets.png')


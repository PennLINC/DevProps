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

% initialize feature weight matrix for SSPs and vectors
SSP_featWeights=zeros(40,174095);
Vec_featWeights=zeros(40,40960);

% reconstruct into original feature space, by adding the absolute feature weight of each PC in accordance with its PC loading


% for each principal component
for P=1:40
	% extract PC weights from this PLS component
	SSP_PCws=SSPWeights(P,:,:);
	% average over all iterations
	av_SSP_PCws=mean(mean((SSP_PCws)));
	% same for vectors
	Vec_PCws=VecWeights(P,:,:);
	% average over all iterations
	av_Vec_PCws=mean(mean((Vec_PCws)));
	% multiply by PC weights to get into feature space
	% pc weights in abs.
	SSP_featWeights(P,:)=(av_SSP_PCws*abs(SSP_PCA.coeff(:,P)));
	Vec_featWeights(P,:)=(av_Vec_PCws*abs(Vec_PCA.coeff(:,P)));
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
        Kend_f=Kstart_f+(5120*2);
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

% print TD BU horz and vert

% INCOMPLETE past this point, need to check 174095 vs 174096 corresp.

% extract out each net and vis 
for k=1:17
        % get indices belonging to this network
        Inds=Kind_f{k};
        % extract the left face for this network
        NetIndsL=Inds(1:5120);
        FaceVecL=coeff(NetIndsL,compNum);
        % extract the right face for this network
        NetIndsR=Inds(5121:10240);
        FaceVecR=coeff(NetIndsR,compNum);
        % print mean coefficient loadings onto this network
        mean(abs(FaceVecR))
        % create png on surface
        Fn=['~/results/PWs/Comp' num2str(compNum) '_Net' num2str(k) '.png'];
        Vis_FaceVec(FaceVecL,FaceVecR,Fn)
end



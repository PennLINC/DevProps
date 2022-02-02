TopographyMat=load('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjTopog.mat');
TopographyMat=TopographyMat.FaceTopogMat;

% add paths
addpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox')

% note transpose
[coeff,score,latent,tsquared,explained,mu]=pca(TopographyMat');

% print out var explained by some components
explained(1:40)

%%% BETTA SAVE DEM PCS BOI

SSP_PC_struct=struct;
SSP_PC_struct.coeff=coeff;
SSP_PC_struct.score=score;
SSP_PC_struct.latent=latent;
SSP_PC_struct.tsquared=tsquared;
SSP_PC_struct.explained=explained;
SSP_PC_struct.mu=mu;

save('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjTopogPCA.mat','SSP_PC_struct')

%%%%% map coefficients back onto the 17 networks

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


% for this component
compNum=1

% print out a PBP for each networks
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

VecMat=load('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjVecComps.mat');
VecMat=VecMat.FaceVecCompMat;

% add paths
addpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox')

% note transpose
[coeff,score,latent,tsquared,explained,mu]=pca(VecMat(1:40960,:)');

% print out var explained by some components
explained(1:40)

%%% BETTA SAVE DEM PCS BOI

Vecs_PC_struct=struct;
Vecs_PC_struct.coeff=coeff;
Vecs_PC_struct.score=score;
Vecs_PC_struct.latent=latent;
Vecs_PC_struct.tsquared=tsquared;
Vecs_PC_struct.explained=explained;
Vecs_PC_struct.mu=mu;

save('/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjVecsPCA.mat','Vecs_PC_struct')

%%%%% map coefficients back onto the sep. components

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

compNum=1

% print out a PBP for each directional component
% extract the left face for this network
FaceVecL=coeff(BULH,compNum);
% extract the right face for this network
FaceVecR=coeff(BURH,compNum);
% print mean coefficient loadings onto this network
mean(abs(FaceVecR))
% create png on surface
Fn=['~/results/PWs/Comp' num2str(compNum) '_BUH.png'];
Vis_FaceVec(FaceVecL,FaceVecR,Fn)

% extract the left face for this network
FaceVecL=coeff(BULV,compNum);
% extract the right face for this network
FaceVecR=coeff(BURV,compNum);
% print mean coefficient loadings onto this network
mean(abs(FaceVecR))
% create png on surface
Fn=['~/results/PWs/Comp' num2str(compNum) '_BUV.png'];
Vis_FaceVec(FaceVecL,FaceVecR,Fn)

% extract the left face for this network
FaceVecL=coeff(TDLH,compNum);
% extract the right face for this network
FaceVecR=coeff(TDRH,compNum);
% print mean coefficient loadings onto this network
mean(abs(FaceVecR))
% create png on surface
Fn=['~/results/PWs/Comp' num2str(compNum) '_TDH.png'];
Vis_FaceVec(FaceVecL,FaceVecR,Fn)

% extract the left face for this network
FaceVecL=coeff(TDLV,compNum);
% extract the right face for this network
FaceVecR=coeff(TDRV,compNum);
% print mean coefficient loadings onto this network
mean(abs(FaceVecR))
% create png on surface
Fn=['~/results/PWs/Comp' num2str(compNum) '_TDV.png'];
Vis_FaceVec(FaceVecL,FaceVecR,Fn)

% extract the left face for this network
FaceVecL=coeff(PropBUL,compNum);
% extract the right face for this network
FaceVecR=coeff(PropBUR,compNum);
% print mean coefficient loadings onto this network
mean(abs(FaceVecR))
% create png on surface
Fn=['~/results/PWs/Comp' num2str(compNum) '_PropBU.png'];
Vis_FaceVec(FaceVecL,FaceVecR,Fn)

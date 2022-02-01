% addpaths

%%%%%%%%%%%%%%%%%%%%%%% Aggregate Subj Loadings
% get subj list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);

% initialize matrix (2 hemis, 17 networks bb)
numFaces=(5120*2)*17;
numSubjs=height(subjs);
InitMat=zeros(numFaces,numSubjs);

% for each subj
for s=1:numSubjs
	subjcell=table2array(subjs(s,1));
        subj=subjcell{:};
	% for each network
	for k=1:17
		subjFPL='/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/resamp_3k_network_L' num2str(k) '.func.gii';
		subjFPR='/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/resamp_3k_network_R$' num2str(k) '.func.gii';
		% load in mat
		% plop in InitMat
	end
end
% end
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% convert to faces
% convert vertex ts to face ts
sizeTsl=size(TRs_l);
lengthTSl=sizeTsl(2);
sizeTsr=size(TRs_r);
lengthTSr=sizeTsr(2);
% initialize face TS : # of faces and # of TRs
faceTS_l=zeros(length(faces_l),sizeTsl(2));
faceTS_r=zeros(length(faces_r),sizeTsr(2));
% there is probably a more efficient way to do this than loop over every TR
for t=1:lengthTSl
        ts_this_tr_l=TRs_l(:,t);
        faceTS_l(:,t)=sum(ts_this_tr_l(faces_l),2)./3;
        ts_this_tr_r=TRs_r(:,t);
        faceTS_r(:,t)=sum(ts_this_tr_r(faces_r),2)./3;
end
%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Run PCA across subjs
fastpca
%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%% print out some PCs
% subj PC values
% visualization of PCA loadings
%%%%%%%%%%%%%%%%%%%%%%%%

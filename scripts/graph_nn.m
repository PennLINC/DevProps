function louv(subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Load in FSC, find optimal community structure, calculate within comm connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Need the brain connectivity toolbox for this one
addpath('~/multiscale/scripts/derive_parcels/Toolbox/2017_01_15_BCT')
% Make gamma distribution for louvain parameter sweep
rng('default');
gammas=.2:.2:6;
% load in subj data
fn=['~/results/PWs/Proced/' subj '/' subj '_FSC.mat'];
sdata=load(fn);
% pull out left and right W matrices
WL=sdata.adjMats.L;
WR=sdata.adjMats.R;
% normalize to 1 is max
WL=WL./max(WL);
WR=WR./max(WR);
% remove under 90% of max
thresh=prctile(WL,90,1);
for F=1:length(WL)
	RowF=WL(F,:);
	surviVec=RowF>thresh(F);
	WL(F,:)=WL(F,:).*surviVec;
end
% right
thresh=prctile(WR,90,1);
for F=1:length(WR)
        RowF=WR(F,:);
        surviVec=RowF>thresh(F);
        WR(F,:)=WR(F,:).*surviVec;
end

%%%%%%%%%%%%%%%%%%%%%% initialize output structures
% community affiliation matrices
AffilMat_L=zeros(30,length(WL));
AffilMat_R=zeros(30,length(WR));
% within_fsc vector for continuous, non community ID plotting
withinVec_L=zeros(1,length(WL));
withinVec_R=zeros(1,length(WR));

%%%%%%%%%%%%%%%%%%%%% for left
disp('left')
%%% louv a few times
for lou=1:30;
	% get community partition
	[M,Q]=community_louvain(WL,1);
	AffilMat_L(lou,:)=M;
end
% get agreement matrix
agMatL=agreement(AffilMat_L');
% replace affilmat with louvain on agreement
for lou=1:30;
        % get community partition
        [M,Q]=community_louvain(agMatL,1);
        AffilMat_L(lou,:)=M;
end
% get max ARI
adj_rand_vec=zeros(1,30);
adj_rand_mat=zeros(30,30);
for arIter=1:30
        arIter
        % sep. loop for comparison to all other 20 solutions
        for r=1:30
                adj_rand_mat(arIter,r)=ri(AffilMat_L(arIter,:),AffilMat_L(r,:),'adjusted');
        end
        adj_rand_vec(arIter)=nanmean(adj_rand_mat(arIter,:));
end
% choose winning solution
[maximum,ind]=max(adj_rand_vec);
winL=AffilMat_L(ind,:);
% calc within module FSC for each face (in winning solution)
for F=1:length(WL)
	community=winL(F);
	withinVec_L(F)=mean(WL(F,winL==community));
end
%%%%%%%%%%%%%%%%%%%%% for right
disp('right')
%%% louv a few times
for lou=1:30;
        % get community partition
        M=community_louvain(WR);
        AffilMat_R(lou,:)=M;
end

% get agreement matrix
agMatR=agreement(AffilMat_R');
% replace affilmat with louvain on agreement
for lou=1:30;
        % get community partition
        [M,Q]=community_louvain(agMatR,1);
        Q
        AffilMat_R(lou,:)=M;
end
% get max ARI
adj_rand_vec=zeros(1,30);
adj_rand_mat=zeros(30,30);
for arIter=1:30
        arIter
        % sep. loop for comparison to all other 20 solutions
        for r=1:30
                adj_rand_mat(arIter,r)=ri(AffilMat_R(arIter,:),AffilMat_R(r,:),'adjusted');
        end
        adj_rand_vec(arIter)=nanmean(adj_rand_mat(arIter,:));
end
% choose winning solution
[maximum,ind]=max(adj_rand_vec);
winR=AffilMat_R(ind,:);
% calc within module FSC for each face (in winning solution)
for F=1:length(WR)
        community=winR(F);
        withinVec_R(F)=mean(WR(F,winR==community));
end

% save out

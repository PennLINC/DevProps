addpath('~/chead_refugees/BCT/2017_01_15_BCT')
addpath('~/net_neuroScripts')
% for schaef parc iteration
cd('~/chead_refugees/schaef400_694_rest')
d=dir('*.txt');
subjs_lm={d.name};
schaef=1:length(subjs_lm);
for i=1:length(subjs_lm);
    dude=strsplit(char(subjs_lm(i)),'_');
    schaef(i)=str2double(dude(1));
end

% for laus
%cd('/data/jux/BBL/studies/pnc/pncDataFreeze20170905/n1601_dataFreeze/neuroimaging/rest/restNetwork_Lausanne250/Lausanne250Networks')
%d=dir('*.txt');

laus={d.name};

% just to get example strings of interest
a=strsplit(char(laus(5)),'_');

% to retain only laus parcs that are retained in schaef parcs
for i=1:length(subjs_lm);
    b=strsplit(char(subjs_lm(i)),'_');
    C={char(b(1)),char(a(2)),char(a(3))};
    subjs_lm(i)=cellstr(strjoin(C,'_'));
end

% retain only those in schaef atlases via reference above restriction of subjects to just those founds in schaef
% assumedly those with low motion
subjs=subjs_lm;

% addition to catch those missed
subjs=fileread('~/not_ran_sch400rs.txt');
subjs=strsplit(subjs);
subjs=subjs(1:150);

% Make gamma distribution
rng('default');
%distr=makedist('normal');
%distr.mu=1;
%gammas=random(distr,100,1);
gammas=-0:.05:6;

% Set params for BCT
B='negative_asym';

% Make output tables for populating later
modmat=cell((length(subjs)),(length(gammas)+1));
avgrandmat=cell((length(subjs)),(length(gammas)+1));
maxrqrandmat=cell((length(subjs)),(length(gammas)+1));
winmat=cell((length(subjs)),(length(gammas)+1));
bwmat=cell((length(subjs)),(length(gammas)+1));
numcommsmat=cell((length(subjs)),(length(gammas)+1));
affilmat=cell((length(gammas)),(400+1));


% loop through subjs
for x=1:10
%for x=1:length(subjs)
    
	mod=[(subjs(x))];
	%seg=[(subjs(x))];
	avgrand=[(subjs(x))];
	maxrqrand=[(subjs(x))];
	win=[(subjs(x))];
	bw=[(subjs(x))];
	numcomms=[(subjs(x))];
	affilmat=zeros(length(gammas),400);

	W = load(char(subjs(x)));
    
    
    % deal with NaNs
    %nans=length(W(isnan(W)));
    %W(isnan(W))=0;
    
	% loop through gammas
    	for g=1:length(gammas);

	% 12/13 addition - 20 runs for each subj at each gamma, extract highest modularity and most consistent(RAND-evaled). Note if they are the same
		win_gamma_solutions=zeros(400,20);
		win_gamma_qs=zeros(1,20);
		avg_adj_rand_vec=zeros(1,20);	
		% get 20 solutions
		for z=1:20;
	    		[M,Q]=community_louvain(W,gammas(g),[],B);
			win_gamma_solutions(:,z)=M;
			win_gamma_qs(z)=Q;
		end
		% get solution with highest modularity
		maxmod=max(win_gamma_qs);
		maxqind=find(win_gamma_qs==maxmod);
		% just take one in case a few are the same partition
		winner_q=maxqind(1);
		% get solution with most between-solution similarity via adjusted rand index
		for r=1:20;
			adj_rand_vec=zeros(1,20);
			% sep. loop for comparison to all other 20 solutions
			for t=1:20;
			adj_r=ri(win_gamma_solutions(:,r),win_gamma_solutions(:,t),'adjusted');
			adj_rand_vec(t)=adj_r;
			end
			avg_adj_rand_vec(r)=nanmean(adj_rand_vec);
		end
		maxrand=nanmax(avg_adj_rand_vec);
		maxrandind=find(avg_adj_rand_vec==maxrand);
		% just take one in case top solutions are the same
		winner_r=maxrandind(1);

		% how similar are maxrand and maxq?
		q_r_sim=ri(win_gamma_solutions(:,winner_r),win_gamma_solutions(:,winner_q),'adjusted');
		
		%%%% will just take rand winner for this iteration, can rerun and use maxmod instead %%%
		M=win_gamma_solutions(:,winner_r);
			
		% Define Modules and Nodes in network
		unique_S=unique(M);
		numNodes=length(W);

		% Number of communities 
		numComm=length(unique_S);

		% Set diagonal of adjacency matrix to nan
		%W=W + diag(repmat(nan,[numNodes,1]));

		% Define community by community matrix
		comm_comm_mat=zeros(numComm,numComm);

		% Define Within/Between Module connectivity matrix
		comm_wb_mat=zeros(numComm,2);
		wb_vec=zeros(1,2);
		com1 = 1;
		for i=unique_S'
			com2 = 1;
			% Define index for nodes in each community
			comidx = find(M==i);
			not_comidx=find(M~=i);
			for j = unique_S'
				comidx_2= find(M==j);
				% Get mean edge weights of edges connecting each pair of communities
				% Pair-wise Between-module connectivity
				current_edges=W(comidx,comidx_2);
				mean_edgeWeight=nanmean(nanmean(current_edges));
				% Define a community X community matrix for each pair of communities
				comm_comm_mat(com1,com2)=mean_edgeWeight;
				com2= com2 + 1;
			end

			% Within module connectivity
			comm_wb_mat(i,1) = nanmean(nanmean(W(comidx,comidx)));
			% Between module connectivity
			comm_wb_mat(i,2) = nanmean(nanmean(W(comidx,not_comidx)));

			com1 = com1 + 1;

		end

		% Compute the overall average within- and between-module connectivity
		within = logical(bsxfun(@eq,M,M'));
		wb_vec(1) = nanmean(W(within));
		wb_vec(2) = nanmean(W(~within));

		% it done don't work the same for dem neggy weighters
		%within_between_ratio = wb_vec(1) / wb_vec(2);

		% Average Within-Module Connectivity
		Avg_Within_Conn=wb_vec(1);
		% Average Between-Module Connectivity
		Avg_Between_Conn=wb_vec(2);
		% Participation Coefficient
		%[pos,neg]=participation_coef_sign(W,M);	
	    
		% print out iteration info
		info = [subjs(x), g, gammas(g), numComm, mean(avg_adj_rand_vec), q_r_sim]
		Metrics = [char(subjs(x)), ', ' num2str(wb_vec(1)), ', ' num2str(wb_vec(2)), ',' num2str(Q), ',' num2str(q_r_sim)];
		mod=[mod;num2str(Q)];
		%seg=[seg;num2str(within_between_ratio)];
		%clustp=[clustp;num2str(mean(pos))];
		%clustn=[clustn;num2str(mean(neg))];
		win=[win;num2str(wb_vec(1))];
		bw=[bw;num2str(wb_vec(2))];
		numcomms=[numcomms;num2str(numComm)];
		avgrand=[avgrand;num2str(mean(avg_adj_rand_vec))];
		maxrqrand=[maxrqrand;num2str(q_r_sim)];
		% put em in da tables, 1:g+1 adds a spot for subj name
		modmat(x,1:g+1)=mod;
		%segmat(x,1:g+1)=seg;
		avgrandmat(x,1:g+1)=avgrand;
		maxrqrandmat(x,1:g+1)=maxrqrand;
		winmat(x,1:g+1)=win;
		bwmat(x,1:g+1)=bw;
		numcommsmat(x,1:g+1)=numcomms;
	    	affilmat(g,:)=M';
	    
	    end
	    csv_name=strcat('~/affilmats/',subjs(x),'_af.csv');
	    csvwrite(char(csv_name),affilmat);
end
csv_name2=strcat('~/numcomsmats/',x,'_nc.csv');
csvwrite(char(csv_name2),numcommsmat);
csv_name3=strcat('~/randstuff/',x,'_ar.csv');
csvwrite(char(csv_name3),avgrandmat);
csv_name4=strcat('~/randstuff/',x,'_rqr.csv');
csvwrite(char(csv_name4),maxrqrandmat);
csv_name5=strcat('~/connecmats/',x,'_win.csv');
csvwrite(char(csv_name5),winmat);
csv_name6=strcat('~/connecmats/',x,'_bw.csv');
csvwrite(char(csv_name6),bwmat);

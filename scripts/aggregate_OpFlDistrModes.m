% load subjs list
Subjs=readtable('~/PWs/hcpd_subj_list.txt')
% initialize intermed array (to take another set of modes from)
IntermedArray_L=zeros(20480,height(Subjs));
IntermedArray_R=zeros(20480,height(Subjs));

%initialize output array
OutArray_L=zeros(20480,1);
OutArray_R=zeros(20480,1);

% edges to apply in discretize
edges = 0:10:180 ;

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,1));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_AngDistMat.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_AngDistMat.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		%discretize each face's distribution
		for f = 1:20480
			Angbins_L=discretize(AngsL(:,f),edges);
			%get mode of distribution
			IntermedArray_L(f,s)=mode(Angbins_L);
			% and for right
			Angbins_R=discretize(AngsR(:,f),edges);
			IntermedArray_R(f,s)=mode(Angbins_R);
		end	
	end
end

% find where subjects were missing, remove
PopulatedColsL=find(sum(IntermedArray_L)~=0);
PopulatedColsR=find(sum(IntermedArray_R)~=0);
IntermedArray_L=IntermedArray_L(:,PopulatedColsL);
IntermedArray_R=IntermedArray_R(:,PopulatedColsR);

%get modes across subjs
for f = 1:20480
	OutArray_L(f,1)=mode(IntermedArray_L(f,:));
	OutArray_R(f,1)=mode(IntermedArray_R(f,:));	
end

%save
writetable(table(OutArray_L),'/cbica/projects/pinesParcels/results/PWs/ModeModes_L.csv');
writetable(table(OutArray_R),'/cbica/projects/pinesParcels/results/PWs/ModeModes_R.csv');


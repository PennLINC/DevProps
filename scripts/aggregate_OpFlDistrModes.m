% load subjs list
Subjs=readtable('~/PWs/hcpd_subj_list.txt')
% initialize intermed array (to take another set of modes from)
IntermedArray=zeros(20484,height(Subjs));
%initialize output array
OutArray=zeros(20480,1);

% edges to apply in discretize
edges = -180:15:180 ;

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,1));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_AngDistMat.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_AngDistMat.mat']);
		AngsL=Angs.AngDist.Left;
		%discretize each face's distribution
		for f = 1:20480
			Angbins=discretize(AngsL(:,f),edges);
			%get mode of distribution
			IntermedArray(f,s)=mode(Angbins);
		end	
	end
end

% find where subjects were missing, remove
PopulatedCols=find(sum(IntermedArray)~=0);
IntermedArray=IntermedArray(:,PopulatedCols);

%get modes across subjs
for f = 1:20480
OutArray(f,1)=mode(IntermedArray(f,:));
end

%save
writetable(table(OutArray),'/cbica/projects/pinesParcels/results/PWs/ModeModes_L.csv');

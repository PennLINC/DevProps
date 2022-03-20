% load subjs list
Subjs=readtable('~/PWs/hcpd_subj_list.txt')
% initialize intermed array (to take another set of modes from)
IntermedArray_L=zeros((5120*4),height(Subjs));
IntermedArray_R=zeros((5120*4),height(Subjs));

%initialize output array
OutArray_L=zeros((5120*4),1);
OutArray_R=zeros((5120*4),1);

% edges to apply in discretize
edges = 0:10:180;

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,1));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_curvAngDistMat.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_curvAngDistMat.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		%discretize each face's distribution
		for f = 1:(5120*4)
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
for f = 1:(5120*4)
	OutArray_L(f,1)=mode(IntermedArray_L(f,:));
	OutArray_R(f,1)=mode(IntermedArray_R(f,:));	
end

%save
writetable(table(OutArray_L),'/cbica/projects/pinesParcels/results/PWs/ModeModes5_curv_L.csv');
writetable(table(OutArray_R),'/cbica/projects/pinesParcels/results/PWs/ModeModes5_curv_R.csv');

% for running on single subjs
subj='';
fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
sOutArray_L=zeros((5120*4),1);
sOutArray_R=zeros((5120*4),1);
% load in subj's distr
Angs=load([fp '/' subj '_AngDistMat5.mat']);
AngsL=Angs.AngDist.gLeft;
AngsR=Angs.AngDist.gRight;
faceModesL_Prom=zeros(1,(5120*4));
faceModesR_Prom=zeros(1,(5120*4));
%discretize each face's distribution
for f = 1:(5120*4)
	%%%% find mode and prominence
	binnedFaceDirs=histcounts(AngsL(:,f),edges);	
	[M,I]=max(binnedFaceDirs);
	sOutArray_L(f)=I;
	%%% find prominence from non-adjacent modes
	% make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
	tmpBFDvec=zeros(1,22);
	tmpBFDvec(3:20)=binnedFaceDirs;
	% I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
	tmpBFDvec((I):(I+4))=[];
	% get updated binnedFaceDirs: those non-adjacent
	[M2,I2]=max(tmpBFDvec);
	% get relative prominence of secondary peak
	faceModesL_Prom(f)=1-(M2/M);
	% and for right
        binnedFaceDirs=histcounts(AngsR(:,f),edges);
        [M,I]=max(binnedFaceDirs);
        sOutArray_R(f)=I;
        %%% find prominence from non-adjacent modes
        % make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
        tmpBFDvec=zeros(1,22);
        tmpBFDvec(3:20)=binnedFaceDirs;
        % I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
        tmpBFDvec((I):(I+4))=[];
        % get updated binnedFaceDirs: those non-adjacent
        [M2,I2]=max(tmpBFDvec);
        % get relative prominence of secondary peak
        faceModesR_Prom(f)=1-(M2/M);
end

Vis_FaceVec_modes5(sOutArray_L,sOutArray_R,'fn2_5.png',faceModesL_Prom,faceModesR_Prom)

% and just include BUprop plotting so its all in the same place for single subj instances

outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs5.mat'];
outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs5.mat'];
% load in L
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,1));
csv_R=table(matR(:,1));

Vis_FaceVec_Bup5(table2array(csv_L),table2array(csv_R),'fnbup_5.png')


% load subjs list
Subjs=readtable('~/PWs/rs_subs.csv')
% initialize intermed array (to take another set of modes from)
IntermedArray_L=zeros(5120,height(Subjs));
IntermedArray_R=zeros(5120,height(Subjs));

%initialize output array
OutArray_L=zeros(5120,1);
OutArray_R=zeros(5120,1);
ModesL_Prom=zeros(5120,height(Subjs));
ModesR_Prom=zeros(5120,height(Subjs));
GroPromArray_L=zeros(5120,1);
GroPromtArray_R=zeros(5120,1);

% edges to discretize on 
edges=0:10:180;

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_curvAngDistMat4.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_curvAngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		%discretize each face's distribution
		for f = 1:5120
			binnedFaceDirs=histcounts(AngsL(:,f),edges);
        		[M,I]=max(binnedFaceDirs);
        		IntermedArray_L(f,s)=I;
        		%%% find prominence from non-adjacent modes
        		% make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
        		tmpBFDvec=zeros(1,22);
        		tmpBFDvec(3:20)=binnedFaceDirs;
        		% I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
        		tmpBFDvec((I):(I+4))=[];
        		% get updated binnedFaceDirs: those non-adjacent
        		[M2,I2]=max(tmpBFDvec);
       			 % get relative prominence of secondary peak
        		ModesL_Prom(f,s)=1-(M2/M);
        		% and for right
        		binnedFaceDirs=histcounts(AngsR(:,f),edges);
        		[M,I]=max(binnedFaceDirs);
        		IntermedArray_R(f,s)=I;
        		%%% find prominence from non-adjacent modes
        		% make a temp. binnedFaceDirs vector with padding on sides: will allow for uniform neighbor removal
        		tmpBFDvec=zeros(1,22);
        		tmpBFDvec(3:20)=binnedFaceDirs;
        		% I and I+4 instead of I-2 I+2 because we added 2 0s to the start of the vector (-2 -> 0,+2 -> 4)
        		tmpBFDvec((I):(I+4))=[];
        		% get updated binnedFaceDirs: those non-adjacent
        		[M2,I2]=max(tmpBFDvec);
        		% get relative prominence of secondary peak
        		ModesR_Prom(f,s)=1-(M2/M);
		end	
	end
end

% find where subjects were missing, remove
PopulatedColsL=find(sum(IntermedArray_L)~=0);
PopulatedColsR=find(sum(IntermedArray_R)~=0);
IntermedArray_L=IntermedArray_L(:,PopulatedColsL);
IntermedArray_R=IntermedArray_R(:,PopulatedColsR);
ModesL_Prom=ModesL_Prom(:,PopulatedColsL);
ModesR_Prom=ModesR_Prom(:,PopulatedColsR);

%get modes across subjs, avg. prominence
for f = 1:5120
	OutArray_L(f,1)=mode(IntermedArray_L(f,:));
	OutArray_R(f,1)=mode(IntermedArray_R(f,:));
	GroPromArray_L(f,1)=mean(ModesL_Prom(f,:));
	GroPromArray_R(f,1)=mean(ModesR_Prom(f,:));	
end

%save
writetable(table(OutArray_L),'/cbica/projects/pinesParcels/results/PWs/curvModeModes4_L.csv');
writetable(table(OutArray_R),'/cbica/projects/pinesParcels/results/PWs/curvModeModes4_R.csv');

Vis_FaceVec_modes(OutArray_L,OutArray_R,'GroupCGModes',GroPromArray_L,GroPromArray_R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for running on single subjs
subj='';
fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
sOutArray_L=zeros(5120,1);
sOutArray_R=zeros(5120,1);
% load in subj's distr
Angs=load([fp '/' subj '_curvAngDistMat4.mat']);
AngsL=Angs.AngDist.gLeft;
AngsR=Angs.AngDist.gRight;
faceModesL_Prom=zeros(1,5120);
faceModesR_Prom=zeros(1,5120);
%discretize each face's distribution
for f = 1:5120
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

Vis_FaceVec_modes(sOutArray_L,sOutArray_R,'fn2_4.png',faceModesL_Prom,faceModesR_Prom)

% and just include BUprop plotting so its all in the same place for single subj instances

outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_curv4.mat'];
outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_curv4.mat'];
% load in L
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,1));
csv_R=table(matR(:,1));

Vis_FaceVec_Bup(table2array(csv_L),table2array(csv_R),'fnbup_4.png')


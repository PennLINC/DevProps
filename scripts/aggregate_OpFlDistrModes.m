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
% and a top-down proportion vector for left and right
TDP_mL=zeros(5120,height(Subjs));
TDP_mR=zeros(5120,height(Subjs));
TDP_vL=zeros(5120,1);
TDP_vR=zeros(5120,1);
TDP_sL=zeros(5120,1);
TDP_sR=zeros(5120,1);

% edges to discretize on 
edges=0:10:180;

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_AngDistMat4.mat']);
		s
		% load in subj's distr
		Angs=load([fp '/' Subj{:} '_AngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% and BUTD vecs
		outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:} '/' Subj{:} '_BUTD_L_resultantVecs.mat'];
		outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:} '/' Subj{:} '_BUTD_R_resultantVecs.mat'];
		% load in L
		BUTD_L=load(outFP_L);
		BUTD_L=cell2mat(BUTD_L.OutDf_L);
		% load in R
		BUTD_R=load(outFP_R);
		BUTD_R=cell2mat(BUTD_R.OutDf_R);
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
			% insert prop. TDm note contingencies. BUTD file is already masked
			if f < 4590;
				TDP_mL(f,s)=1-BUTD_L(f,1);
			end
			if f < 4596
				TDP_mR(f,s)=1-BUTD_R(f,1);
			end
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
	% same contingencies as above
	if f < 4590;
		TDP_vL(f,1)=mean(TDP_mL(f,:));
		TDP_sL(f,1)=std(TDP_mL(f,:));
	end
	if f < 4596
		TDP_vR(f,1)=mean(TDP_mR(f,:));	
		TDP_sR(f,1)=std(TDP_mR(f,:));
	end
end

% read magnitude length vecs for plot masking
magL=readtable('~/data/vHL.csv');
magR=readtable('~/data/vHR.csv');

%save
writetable(table(OutArray_L),'/cbica/projects/pinesParcels/results/PWs/ModeModes4_L.csv');
writetable(table(OutArray_R),'/cbica/projects/pinesParcels/results/PWs/ModeModes4_R.csv');

Vis_FaceVec_modes(OutArray_L,OutArray_R,'GroupPGGModes',GroPromArray_L,GroPromArray_R)

% filling in unmasked faces with unmasked BU Prop averages, which are the first 45xx
Vis_FaceVec_Bup(TDP_vL(1:4589),TDP_vR(1:4595),'GroupTDProp.png')
Vis_FaceVec_Bup_desat(TDP_vL(1:4589),TDP_vR(1:4595),'GroupTDProp_desat.png',table2array(magL),table2array(magR))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for running on single subjs
subj='';
fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
sOutArray_L=zeros(5120,1);
sOutArray_R=zeros(5120,1);
% load in subj's distr
Angs=load([fp '/' subj '_AngDistMat4.mat']);
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

outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
% load in L
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(1-matL(:,1));
csv_R=table(1-matR(:,1));

Vis_FaceVec_Bup(table2array(csv_L),table2array(csv_R),'ExampleSubjTDprop.png')


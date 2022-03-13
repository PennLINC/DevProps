% add paths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% read in subj list
Subjs=readtable('~/PWs/rs_subs.csv');
%initialize output array (one row for subj name, one for dmn segreg
DMNseg=cell(height(Subjs),2);
% set parent FP for shorter FPs in loop
parentfp = '/cbica/projects/hcpd/data/motMasked_contSegs/';
sspParentfp='/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/';
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';

% get medial wall mask
% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);
nomw_mask_l=setdiff(1:2562,mw_v_l);
nomw_mask_r=setdiff(1:2562,mw_v_r);

% for each subj
for s=1:height(Subjs)
	subj=table2array(Subjs{s,2});
	% load in TRs_l and TRs_r
	TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_3k.mgh'];
	TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_3k.mgh'];
	% filepaths to files
	TRs_lf=MRIread(TRs_lfp);
	TRs_rf=MRIread(TRs_rfp);
	% files to data
	TRs_l=squeeze(TRs_lf.vol);
	TRs_r=squeeze(TRs_rf.vol);
	% ts mask - 10 extra MW verts induces from .mgh transformation
	MW_ind_l=find(std(TRs_l')==0);
	MW_ind_r=find(std(TRs_r')==0);
	% combine medial wall label with mw .mgh
	mw_v_mask_l=union(mw_v_l,MW_ind_l);
	mw_v_mask_r=union(mw_v_r,MW_ind_r);
	% get non mw indices
	nomw_mask_l=setdiff(1:2562,mw_v_mask_l);
	nomw_mask_r=setdiff(1:2562,mw_v_mask_r);
	% mergedTS, with mw masked out
	mTS=vertcat(TRs_l(nomw_mask_l,:),TRs_r(nomw_mask_r,:));
	% initialize network partition-holder
	PartFrame=zeros((length(nomw_mask_l)+length(nomw_mask_r)),17);
	% load in network partition
	for k=1:17
		subjFPL=[sspParentfp subj '/resamp_3k_network_L' num2str(k) '.func.gii'];
	        subjFPR=[sspParentfp subj '/resamp_3k_network_R' num2str(k) '.func.gii'];
		V_Topog_L=gifti(subjFPL);
        	V_Topog_R=gifti(subjFPR);
		% extract loadings from file, plant into vertices by networks matrix for this subject
       		PartFrame(1:length(nomw_mask_l),k)=V_Topog_L.cdata(nomw_mask_l,1);
        	PartFrame((length(nomw_mask_l)+1):(length(nomw_mask_l)+length(nomw_mask_r)),k)=V_Topog_R.cdata(nomw_mask_r,1);
	end
	% take max of each vertex to get hard parcel
	[ ~ , hardParcel]=max(PartFrame,[],2);	
	% initialize FCmat
	FC=zeros(17,17);
	% fc matrix at fs4 resolution
	conmat=corrcoef(mTS');
	for k=1:17
		% get vertices in this network
		VertsInNet=find(hardParcel==k);
		% index out within network connectivity
		curNetMat=conmat(VertsInNet,VertsInNet);
		% get within-network connectivity: using TRIU to get a mask of the upper triangle and avoid the diagonal/redundancy
		wincon=mean(curNetMat(find(~triu(ones(size(curNetMat))))));
		FC(k,k)=wincon;
		% make vector for all values except for current K (N) to loop through
		Kvec=1:17;
		NotKvec=Kvec(Kvec~=k); 
		% one fewer networks than number of networks, as we are looping over OTHER networks (other than k)
		for b=1:16
				curOtherNet=NotKvec(b);
				% index vertices not in up-one-level-network-N loop
				NotKind=find(hardParcel==curOtherNet);
				bwMat=conmat(VertsInNet,NotKind);
				bwcon=mean(mean(bwMat));
				% plop it into 17x17 matrix
				FC(k,curOtherNet)=bwcon;
		end
	end
	% save out FC mat as .mat
	outFN=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_FCmat.mat'];
	save(outFN,'FC')
	% network 11 is the transmodal dmn
	notDMvec=Kvec(Kvec~=11);
	Dseg=mean(FC(11,notDMvec))
	DMNseg(s,1)=cellstr(subj);
	DMNseg(s,2)=num2cell(Dseg);	
end
% save out aggregate dmn 
writetable(table(DMNseg),'~/Results/DMNseg.csv')

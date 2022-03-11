% add paths
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% read in subj list
Subjs=readtable('~/PWs/rs_subs.csv')
%initialize output array (one row for subj name, one for dmn segreg
DMNseg=cell(height(Subjs),2);
% set parent FP for shorter FPs in loop
parentfp = '/cbica/projects/hcpd/data/motMasked_contSegs/';

% for each subj
for s=1:height(Subjs)
	subj=table2array(Subjs{s,2});
	% read in scrubbed time series
	TSfp=['/cbica/projects/hcpd/data/motMasked_contSegs/' subj '/' subj '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii'];
	TS=read_cifti(TSfP);
	% isolate cortex
	TSc=TS.cdata(1:59412,:);
	% read in network partition
	netPart=['/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat']; 
	networks=load(netPart);
	% conver to hard parcels
	softParcel=networks.V{:};
	[ ~ , hardParcel]=max(softParcel,[],2); 
	% initialize FCmat
	FC=zeros(17,17);
	% the return of bigass connectivity matrix: 59k edition
	ba_conmat=corrcoef(TSc);
	for k=1:17
		% get vertices in this network
		VertsInNet=find(hardParcel==k);
		% index out within network connectivity
		curNetMat=ba_conmat(VertsInNet,VertsInNet);
		% get within-network connectivity: using TRIU to get a mask of the upper triangle and avoid the diagonal/redundancy
		wincon=mean(curNetMat(find(~triu(ones(size(curNetMat))))));
		FC(k,k)=wincon;
		% get fc
		
	% save out FC mat as .mat

	% save out DMN seg 
end

% save out aggregate dmn 

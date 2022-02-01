addpath(genpath('~/dropbox/cifti-matlab/'));
% use 600 TRs plus subj list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G600TRs.txt','ReadVariableNames',false);
% read in hard parcel to replace Cort data with
HP=cifti_read('/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients.dscalar.nii');

% loop over each subj
for s=1:height(subjs)
	tic
	subjcell=table2array(subjs(s,1));
	subj=subjcell{:}
	Subjfp=['/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj];
	Loading_Mat = load([Subjfp '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat']);
	Loading_Mat = cell2mat(Loading_Mat.V);
	% for all 17 networks
	for i = 1:17
		HP.cdata(1:59412)=Loading_Mat(:,i);
	  	outputfile=['/cbica/projects/hcpd/data/SSPs/SingleParcel_1by1_fslr/' subj '/AtlasLoading_Network_' num2str(i) '.dscalar.nii'];
	  	cifti_write(HP,outputfile)
	end
	toc
end

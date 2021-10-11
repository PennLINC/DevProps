% add matlab path for used functions
addpath(genpath('/cbica/projects/abcdfnets/scripts/code_nmf_cifti/tool_folder'));
% Load in Subj list
fid = fopen('~/clean100.txt'); 
subjs=textscan(fid, '%s', 'HeaderLines', 0);
subjs=subjs{:};
% initialize output array, 17 nets by 100 subjs (17 -> 34 to record both pgs)
OutArray=zeros(length(subjs),34);
% For S in subjs
for s = 1:length(subjs)
	% load in gradients
	ProjectFolder =['/cbica/projects/abcdfnets/results/wave_output/' char(subjs(s)) '/'];
	PG=read_cifti([ProjectFolder char(subjs(s)) '_PG_LR_32k.dscalar.nii']);
	% extract unimodal-transmodal gradient
	grad1=PG.cdata(:,1);
	% and vismotor
	grad2=PG.cdata(:,2);
	% load in networks
	parcelLoc=strjoin(['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' string(subjs(s)) '/' string(subjs(s)) '_Parcel.dscalar.nii'],'');
	parcel=read_cifti(parcelLoc);
	Nets=parcel.cdata;
	% initialize hierarchy vectors, can calc. grad distance in R
	PG1vec=zeros(17,1);
	PG2vec=zeros(17,1);
	% for each network
	for N=1:17
		NetInds=find(Nets==N);
		PGVals=grad1(NetInds);	
		PG2Vals=grad2(NetInds);
		PG1vec(N)=mean(PGVals);
		PG2vec(N)=mean(PG2Vals);
	end
	OutArray(s,1:17)=PG1vec;
	OutArray(s,18:34)=PG2vec;
end
% save out aggregate array
writetable(table(OutArray),'~/results/NetsInGrads.csv')


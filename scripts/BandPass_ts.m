function BandPass_ts(subj)

% deep gmless cifti for template
sname=char(subj);
parentfp=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
rsfp=[parentfp sname '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii'];
sstfp=[parentfp sname '_ses-baselineYear1Arm1_task-SST_p2mm_masked.dtseries.nii'];
nbackfp=[parentfp sname '_ses-baselineYear1Arm1_task-nback_p2mm_masked.dtseries.nii'];
midfp=[parentfp sname '_ses-baselineYear1Arm1_task-MID_p2mm_masked.dtseries.nii'];

% construct bandpass filter
tr = .8; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
hp_thresh = .01; % lower bound
lp_thresh = .1; % higher bound
[b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));

% set up 
missingDir=['/cbica/projects/abcdfnets/results/MissingDataReports/' sname];
mkdir(missingDir);
missingFile=[missingDir '/MissingData.txt'];

if ~exist(sstfp,'file')
disp('missing sst')
system(['echo sst_missing >> ' missingFile]); 
end

if ~exist(rsfp,'file')
disp('missing rs')
system(['echo rs_missing >> ' missingFile]);
end

if ~exist(nbackfp,'file')
disp('missing nback')
system(['echo nback_missing >> ' missingFile]);
end

if ~exist(midfp,'file')
disp('missing mid')
system(['echo mid_missing >> ' missingFile]);
end


% isolate the masked grayordinatewise time series
% and stack them onto an init oordinate-wise value
Oords=zeros(91282,1);

%%%%%%% RESTING STATE
if exist(rsfp,'file')
rs=read_cifti(rsfp);
rsts=rs.cdata;
% bandpass the data
for v = 1:91282
rsts(v,:)=filtfilt(b,a,double(rsts(v,:)));
end
% reset the cdata matrix to be the clean time series
rs.cdata=rsts;
fpr=[parentfp sname '_p2mm_masked_filtered_rest.dtseries.nii'];
% save out cifti
write_cifti(rs,fpr);
% load in GS to filter
gsfpr=[parentfp sname '_p2mm_masked_rest_GS.csv'];
gs_r=load(gsfpr);
gs_r=filtfilt(b,a,gs_r);
% save out GS
gsfprf=[parentfp sname '_p2mm_masked_filtered_rest_GS.csv'];
writetable(array2table(gs_r),gsfprf,'WriteVariableNames',0);
end
%%%%%%%%%%

%%%%%%% STOP-START
if exist(sstfp,'file')
sst=read_cifti(sstfp);
sstts=sst.cdata;
% bandpass the data
for v = 1:91282
sstts(v,:)=filtfilt(b,a,double(sstts(v,:)));
end
sst.cdata=sstts;
% save out cifti
fps=[parentfp sname '_p2mm_masked_filtered_SST.dtseries.nii'];
write_cifti(sst,fps);
% load in GS to filter
gsfps=[parentfp sname '_p2mm_masked_SST_GS.csv'];
gs_s=load(gsfps);
gs_s=filtfilt(b,a,gs_s);
% save out GS
gsfpsf=[parentfp sname '_p2mm_masked_filtered_SST_GS.csv'];
writetable(array2table(gs_s),gsfpsf,'WriteVariableNames',0);
end
%%%%%%%%%%%

%%%%%%% N-BACK
if exist(nbackfp,'file')
nback=read_cifti(nbackfp);
nbackts=nback.cdata;
% bandpass the data
for v = 1:91282
nbackts(v,:)=filtfilt(b,a,double(nbackts(v,:)));
end
nback.cdata=nbackts;
% save out cifti
fpn=[parentfp sname '_p2mm_masked_filtered_nback.dtseries.nii'];
write_cifti(nback,fpn);
% load in GS to filter
gsfpn=[parentfp sname '_p2mm_masked_nback_GS.csv'];
gs_n=load(gsfpn);
gs_n=filtfilt(b,a,gs_n);
% save out GS
gsfpnf=[parentfp sname '_p2mm_masked_filtered_nback_GS.csv'];
writetable(array2table(gs_n),gsfpnf,'WriteVariableNames',0);
end
%%%%%%%%%%

%%%%%%% EMOTION ID
if exist(midfp,'file')
mid=read_cifti(midfp);
midts=mid.cdata;
% bandpass the data
for v = 1:91282
midts(v,:)=filtfilt(b,a,double(midts(v,:)));
end
mid.cdata=midts;
% save out cifti
fpe=[parentfp sname '_p2mm_masked_filtered_mid.dtseries.nii'];
write_cifti(mid,fpe);
% load in GS to filter
gsfpm=[parentfp sname '_p2mm_masked_MID_GS.csv'];
gs_m=load(gsfpm);
gs_m=filtfilt(b,a,gs_m);
% save out GS
gsfpmf=[parentfp sname '_p2mm_masked_filtered_mid_GS.csv'];
writetable(array2table(gs_m),gsfpmf,'WriteVariableNames',0);
end
%%%%%%%%%%

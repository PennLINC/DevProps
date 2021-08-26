function BandPass_ts(subj)

% deep gmless cifti for template
sname=char(subj);
parentfp=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
rsfp=[parentfp sname '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii'];
sstfp=[parentfp sname '_ses-baselineYear1Arm1_task-SST_p2mm_masked.dtseries.nii'];
nbackfp=[parentfp sname '_ses-baselineYear1Arm1_task-nback_p2mm_masked.dtseries.nii'];
midfp=[parentfp sname '_ses-baselineYear1Arm1_task-MID_p2mm_masked.dtseries.nii'];
tasks=["rest","MID","SST","nback"];
% construct bandpass filter
tr = .8; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
hp_thresh = .008; % lower bound
lp_thresh = .09; % higher bound
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

%%%%%%% RESTING STATE
if exist(rsfp,'file')
	% initialize grayOrd TS
	Oords=zeros(91282,1);
	rs=read_cifti(rsfp);
	rsts=rs.cdata;
	% load in continuous segments mask
	rsCSm_fn=[parentfp sname '_ses-baselineYear1Arm1_task-rest_ValidSegments.txt'];
	rsCSm=dlmread(rsCSm_fn);
	% number of segments
	segShape=size(rsCSm);
	numSegs=segShape(1);
	% for each segment
	for s = 1:numSegs
		% select all vertices from this segment
		% first column is start of clean window, second column is length of clean window
		% -1 because duration of clean window includes first TR indexed
		segment=rsts(:,rsCSm(s,1):(rsCSm(s,1)+rsCSm(s,2)-1));
		% bandpass the data
		for v = 1:91282
			segment(v,:)=filtfilt(b,a,double(segment(v,:)));
		end
		% append bandpasses segment onto grayOrds
		Oords = [Oords segment];
	end
	% remove initialization volume of grayOrds
	Oords(:,1)=[];
	% reset the cdata matrix to be the clean time series
	rs.cdata=Oords;
	TSsize=size(Oords);
	% reset temporal dimension to reflect thresholding of short segments
	rs.diminfo{2}.length=TSsize(2);
	fpr=[parentfp sname '_p2mm_masked_filtered_rest.dtseries.nii'];
	% save out cifti
	write_cifti(rs,fpr);
	% load in GS to filter
	gsfpr=[parentfp sname '_p2mm_masked_rest_GS.csv'];
	gs_r=load(gsfpr);
	% initialize gs
	gs_filt=[];
	% filter each segment sep.
	for s = 1:numSegs
		% -1 because duration of clean window includes first TR indexed
		segment=gs_r(rsCSm(s,1):(rsCSm(s,1)+rsCSm(s,2)-1));
		segSignal=filtfilt(b,a,gs_r);
		gs_filt=[gs_filt segSignal];
	end
	% save out GS
	gsfprf=[parentfp sname '_p2mm_masked_filtered_rest_GS.csv'];
	writetable(array2table(gs_filt),gsfprf,'WriteVariableNames',0);
end
%%%%%%%%%%

%%%%%%% STOP-START
if exist(sstfp,'file')
        % initialize grayOrd TS
        Oords=zeros(91282,1);
	sst=read_cifti(sstfp);
	sstts=sst.cdata;
	% load in continuous segments mask
	sstCSm_fn=[parentfp sname '_ses-baselineYear1Arm1_task-SST_ValidSegments.txt'];
	sstCSm=dlmread(sstCSm_fn);
        % number of segments
        segShape=size(sstCSm);
        numSegs=segShape(1);
        % for each segment
        for s = 1:numSegs
                % select all vertices from this segment
                % first column is start of clean window, second column is length of clean window
                % -1 because duration of clean window includes first TR indexed
                segment=sstts(:,sstCSm(s,1):(sstCSm(s,1)+sstCSm(s,2)-1));
                % bandpass the data
                for v = 1:91282
                        segment(v,:)=filtfilt(b,a,double(segment(v,:)));
                end
                % append bandpasses segment onto grayOrds
                Oords = [Oords segment];
	end
        % remove initialization volume of grayOrds
        Oords(:,1)=[];
        % reset the cdata matrix to be the clean time series
        sst.cdata=Oords;
        TSsize=size(Oords);
        % reset temporal dimension to reflect thresholding of short segments
        sst.diminfo{2}.length=TSsize(2);
	% save out cifti
	fps=[parentfp sname '_p2mm_masked_filtered_SST.dtseries.nii'];
	write_cifti(sst,fps);
	% load in GS to filter
	gsfps=[parentfp sname '_p2mm_masked_SST_GS.csv'];
	gs_s=load(gsfps);
        % initialize gs
        gs_filt=[];
        % filter each segment sep.
        for s = 1:numSegs
                % -1 because duration of clean window includes first TR indexed
                segment=gs_s(sstCSm(s,1):(sstCSm(s,1)+sstCSm(s,2)-1));
                segSignal=filtfilt(b,a,gs_s);
                gs_filt=[gs_filt segSignal];
        end	
	% save out GS
	gsfpsf=[parentfp sname '_p2mm_masked_filtered_SST_GS.csv'];
	writetable(array2table(gs_filt),gsfpsf,'WriteVariableNames',0);
end
%%%%%%%%%%%

%%%%%%% N-BACK
if exist(nbackfp,'file')
        % initialize grayOrd TS
        Oords=zeros(91282,1);
	nback=read_cifti(nbackfp);
	nbackts=nback.cdata;
	% load in continuous segments mask
	nbCSm_fn=[parentfp sname '_ses-baselineYear1Arm1_task-nback_ValidSegments.txt'];
	nbCSm=dlmread(nbCSm_fn);
        % number of segments
        segShape=size(nbCSm);
        numSegs=segShape(1);
        % for each segment
        for s = 1:numSegs
                % select all vertices from this segment
                % first column is start of clean window, second column is length of clean window
                % -1 because duration of clean window includes first TR indexed
                segment=nbackts(:,nbCSm(s,1):(nbCSm(s,1)+nbCSm(s,2)-1));
                % bandpass the data
                for v = 1:91282
                        segment(v,:)=filtfilt(b,a,double(segment(v,:)));
                end
                % append bandpasses segment onto grayOrds
                Oords = [Oords segment];
	end
        % remove initialization volume of grayOrds
        Oords(:,1)=[];
        % reset the cdata matrix to be the clean time series
        nback.cdata=Oords;
        TSsize=size(Oords);
        % reset temporal dimension to reflect thresholding of short segments
        nback.diminfo{2}.length=TSsize(2);
	% save out cifti
	fpn=[parentfp sname '_p2mm_masked_filtered_nback.dtseries.nii'];
	write_cifti(nback,fpn);
	% load in GS to filter
	gsfpn=[parentfp sname '_p2mm_masked_nback_GS.csv'];
	gs_n=load(gsfpn);
        % initialize gs
        gs_filt=[];
        % filter each segment sep.
        for s = 1:numSegs
                % -1 because duration of clean window includes first TR indexed
                segment=gs_n(nbCSm(s,1):(nbCSm(s,1)+nbCSm(s,2)-1));
                segSignal=filtfilt(b,a,gs_n);
                gs_filt=[gs_filt segSignal];
        end
	% save out GS
	gsfpnf=[parentfp sname '_p2mm_masked_filtered_nback_GS.csv'];
	writetable(array2table(gs_filt),gsfpnf,'WriteVariableNames',0);
end
%%%%%%%%%%

%%%%%%% EMOTION ID
if exist(midfp,'file')
	% initialize grayOrd TS
        Oords=zeros(91282,1);
	mid=read_cifti(midfp);
	midts=mid.cdata;
	% load in continuous segments mask
	midCSm_fn=[parentfp sname '_ses-baselineYear1Arm1_task-MID_ValidSegments.txt'];
	midCSm=dlmread(midCSm_fn);
        % number of segments
        segShape=size(midCSm);
        numSegs=segShape(1);
        % for each segment
        for s = 1:numSegs
	        % select all vertices from this segment
                % first column is start of clean window, second column is length of clean window
                % -1 because duration of clean window includes first TR indexed
                segment=midts(:,midCSm(s,1):(midCSm(s,1)+midCSm(s,2)-1));
                % bandpass the data
                for v = 1:91282
                        segment(v,:)=filtfilt(b,a,double(segment(v,:)));
                end
                % append bandpasses segment onto grayOrds
                Oords = [Oords segment];
	end
        % remove initialization volume of grayOrds
        Oords(:,1)=[];
        % reset the cdata matrix to be the clean time series
        mid.cdata=Oords;
        TSsize=size(Oords);
        % reset temporal dimension to reflect thresholding of short segments
        mid.diminfo{2}.length=TSsize(2);
	% save out cifti
	fpe=[parentfp sname '_p2mm_masked_filtered_MID.dtseries.nii'];
	write_cifti(mid,fpe);
	% load in GS to filter
	gsfpm=[parentfp sname '_p2mm_masked_MID_GS.csv'];
	gs_m=load(gsfpm);
        % initialize gs
        gs_filt=[];
        % filter each segment sep.
        for s = 1:numSegs
                % -1 because duration of clean window includes first TR indexed
                segment=gs_m(midCSm(s,1):(midCSm(s,1)+midCSm(s,2)-1));
                segSignal=filtfilt(b,a,gs_m);
                gs_filt=[gs_filt segSignal];
        end	
	% save out GS
	gsfpmf=[parentfp sname '_p2mm_masked_filtered_MID_GS.csv'];
	writetable(array2table(gs_filt),gsfpmf,'WriteVariableNames',0);
end
%%%%%%%%%%

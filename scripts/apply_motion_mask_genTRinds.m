function apply_motion_mask_genTRinds(subj)
% this function motion masks the individual runs (non-filtered) and the concatenated (filtered)
% note that both validsegcell_full and validsegcell_trunc are saved out
% validsegcell_full refers to the TR and span of segment in the outlier FD masked images
% validsegcell truncated refers to the TR and span of segment in the outlierFD masked AND short segment masked images
% minproc and FD leave this script matching truncated, proc TS and GS leave this script matching full
% (full for BandPass, Trunc for deriveWaveProps)
% also, save out unthresh to match basis time series of networks

% load in subj
topleveldir='/cbica/projects/hcpd/data/sub-*'
direc=dir(topleveldir);

% initialize empty vector for average length
TRvecNum=[];

%%% for each "task"
%%% tasks=["rest","MID","SST","nback"];

% filepath for TR # output
ProjectFolder = '/cbica/projects/pinesParcels/results/PWs';
ResultantFolder = [ProjectFolder '/PreProc/'];

task='rest';
sname=subj;

% mask concatenated data
fpParent=['/cbica/projects/hcpd/data/' sname '/ses-V1/files/MNINonLinear/Results/'];
fp=[fpParent 'task-' task '_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii'];

% if file exists, run it
if exist(fp,'file')
	% make output dir
	OutDir=[ResultantFolder sname '/'];
	mkdirCommand=['mkdir ' OutDir];
	system(mkdirCommand);
	% read in cifti
	ts_cif=read_cifti(fp);
	ts=ts_cif.cdata;
	% get size
	ciftisize=size(ts);
	numTRs=ciftisize(2);
	% get cortex indices
	C_ind=1:59412;
	% load in mask
	masfp=['/cbica/projects/hcpd/data/motion/fd/' sname '/ses-V1/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/task-' task '_power_2014_FD_only.mat'];
	if exist(masfp,'file')
		mask=load(masfp);
		% get to FD_thresh of .2 mm, corresponds to threshold 21
		maskp2mm=mask.motion_data{21}.frame_removal;
		TRwise_mask=logical(maskp2mm);
		% length of mask corresponds to number of TRs
		% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
		TRwise_mask=~TRwise_mask;
		% remove TRs with corresp. flag
		masked_trs=ts(:,TRwise_mask);
		% reconfig cifti metadata to reflect new number of TRs
		newciftiSize=size(masked_trs);
		newTRnum=newciftiSize(2);
		% setting continuous frame threshold to 13 TRs in a row
		Threshold=13;			% find changepoints in binary bask
		% find changepoints in binary bask
		d = [true, diff(TRwise_mask') ~= 0];
		% index of changepoints
		dInd=find(d);
		% find difference in indices of changepoints (span of mask/non-mask epochs)
		n = diff([dInd, numTRs]); 
		% find which segments correspond to non-mask
		maskValAtChange=TRwise_mask(dInd);
		ContSegments=n(:,maskValAtChange);
		% create list of starting TR and duration of segments uninterupt. by combined mask
		UTSegSize=size(ContSegments);
		UTSegNum=UTSegSize(2);
		UTSegCell=cell(UTSegNum,2);
		% plant in TR start and duration of clean segments
		for i=1:UTSegNum
			UTSegCell(i,2)=num2cell(ContSegments(i));
		end
		% make 1st column start position in .2mm outlier masked sequence
		% (just the start where prev. segment left off, no masked TRs in gaps b/w)
		UTSegCell(1,1)=num2cell(1);
		for i=2:UTSegNum
			UTSegCell(i,1)=num2cell(UTSegCell{i-1,1}+UTSegCell{i-1,2});
		end
		% check that sum of TRs matches field from dcan/midb mask
		allRetainedSegmentTRLengths=UTSegCell(:,2);
		if (sum([allRetainedSegmentTRLengths{:}])==~mask.motion_data{21}.remaining_frame_count)
			error('dcan/midb remaining combined count does not match internal representation')
		end
		% save remaining_combined_count
		remaining_combined_count=mask.motion_data{21}.remaining_frame_count;
		mask.motion_data{21}.remaining_frame_count
		remaining_cmb_fn=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_remainingTRs.mat'];
		save(remaining_cmb_fn, 'remaining_combined_count');
                % find segments with more continuous TRs than threshold
                OverThreshSegments=find(ContSegments>Threshold);
                % sum remaining segments to get included TRs if this thresh chosen
                RemainingTRs=sum(ContSegments(OverThreshSegments))
                % index of which TR valid segments start at
                ValidTRStarts=dInd(maskValAtChange);
		% index out segments greater than TR thresh from UnThreshSegmentCellstruct
		ValidSegCell=UTSegCell(OverThreshSegments,:);
		% save out version for network basis time series: derived from all low Mot non-outlier volumes
		segmentfnUt=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_ValidSegments_Unthr'];
		writetable(cell2table(UTSegCell),segmentfnUt,'WriteVariableNames',0)
		% save out version for bandpassing - matches TRs of dtseries when it is at that point
		segmentfnFu=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_ValidSegments_Full'];
                writetable(cell2table(ValidSegCell),segmentfnFu,'WriteVariableNames',0)
		% adjust ValidSegCell_Trunc to describe in terms of retained TRs only
		ValidSegCell_Trunc=ValidSegCell;
		ValidSegCell_Trunc(1,1)=num2cell(1);
		for s=2:length(OverThreshSegments)
			prevStart=ValidSegCell_Trunc(s-1,1);
			prevLength=ValidSegCell_Trunc(s-1,2);
			NewStart=prevStart{:}+prevLength{:};
			ValidSegCell_Trunc(s,1)=num2cell(NewStart);
		end	
		% save 2-column df indicating start of valid segments and length: matches dtseries at derive_wave
		segmentfnTr=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_ValidSegments_Trunc'];
		writetable(cell2table(ValidSegCell_Trunc),segmentfnTr,'WriteVariableNames',0)
		% overwite diminfo of cifti
		ts_cif.diminfo{2}.length=newTRnum;
		% overwrite TRs for new file
		ts_cif.cdata=masked_trs;
		% set output filepath
		ofp=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_p2mm_masked.dtseries.nii'];
		% There is no reason this should be a requried step
		ofp=convertStringsToChars(ofp);
		% write out motion masked cifti
		write_cifti(ts_cif,ofp);
end
end

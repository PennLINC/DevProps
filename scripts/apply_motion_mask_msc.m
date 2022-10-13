function apply_motion_mask_msc(subj)
% this function motion masks the concatenated runs 
% note that both validsegcell_full and validsegcell_trunc are saved out: only _Trunc used in this workflow
% validsegcell_full refers to the TR and span of segment in the outlier FD masked images
% validsegcell truncated refers to the TR and span of segment in the outlierFD masked AND short segment masked images
% minproc and FD leave this script matching truncated, proc TS and GS leave this script matching full
% (full for BandPass, Trunc for OpFl)

addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% initialize empty vector for average length
TRvecNum=[];

% filepath for output
ProjectFolder = '/cbica/projects/pinesParcels/results/PWs';
ResultantFolder = ['/cbica/projects/pinesParcels/data/msc/motMasked_contSegs/'];

task='rest';
sname=subj;

% for each run
for r=1:10
	% get strings for this run
	runstrings={'01','02','03','04','05','06','07','08','09','10'};
	runstring=runstrings{r}

	% mask concatenated data
	fpParent=['/cbica/projects/pinesParcels/data/msc/dtseries/' sname '/ses-func' runstring '/files/MNINonLinear/Results/'];
	fp=[fpParent 'task-' task '_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii']

	% if file exists, run it
	if exist(fp,'file')
		% make output dir
		OutDir=[ResultantFolder sname '/'];
		mkdirCommand=['mkdir ' OutDir];
		system(mkdirCommand);
		% read in cifti
		fp
		ts_cif=read_cifti(fp);
		ts=ts_cif.cdata;
		% get size
		ciftisize=size(ts);
		numTRs=ciftisize(2);
		% get cortex indices
		C_ind=1:59412;
		% load in mask
		masfp=['/cbica/projects/pinesParcels/data/msc/' sname '/ses-func' runstring '/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/task-' task '_power_2014_FD_only.mat'];
		if exist(masfp,'file')
			mask=load(masfp);
			% get to FD_thresh of .4 mm, corresponds to threshold 41
			maskp2mm=mask.motion_data{41}.frame_removal;
			TRwise_mask=logical(maskp2mm);
			% length of mask corresponds to number of TRs
			% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
			TRwise_mask=~TRwise_mask;
			% remove TRs with corresp. flag
			masked_trs=ts(:,TRwise_mask);
			% reconfig cifti metadata to reflect new number of TRs
			newciftiSize=size(masked_trs);
			newTRnum=newciftiSize(2);
			% setting continuous frame threshold to 10 TRs in a row
			Threshold=10;			% find changepoints in binary bask
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
			% make 1st column start position in .4mm outlier masked sequence
			% (just the start where prev. segment left off, no masked TRs in gaps b/w)
			UTSegCell(1,1)=num2cell(1);
			for i=2:UTSegNum
				UTSegCell(i,1)=num2cell(UTSegCell{i-1,1}+UTSegCell{i-1,2});
			end
			% check that sum of TRs matches field from dcan/midb mask
			allRetainedSegmentTRLengths=UTSegCell(:,2);
			if (sum([allRetainedSegmentTRLengths{:}])==~mask.motion_data{41}.remaining_frame_count)
				error('dcan/midb remaining combined count does not match internal representation')
			end
			% save remaining_combined_count
			remaining_combined_count=mask.motion_data{41}.remaining_frame_count;
			mask.motion_data{41}.remaining_frame_count
			remaining_cmb_fn=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_remainingTRs_run_' runstring '.mat'];
			save(remaining_cmb_fn, 'remaining_combined_count');
        	        % find segments with more continuous TRs than threshold
        	        OverThreshSegments=find(ContSegments>Threshold);
        	        % sum remaining segments to get included TRs if this thresh chosen
        	        RemainingTRs=sum(ContSegments(OverThreshSegments))
        	        % index of which TR valid segments start at
        	        ValidTRStarts=dInd(maskValAtChange);
			% index out segments greater than TR thresh from UnThreshSegmentCellstruct
			ValidSegCell=UTSegCell(OverThreshSegments,:);
		
			% adjust ValidSegCell_Trunc to describe in terms of retained TRs only
			% in other words, epoch start points should be one "TR" after 
			ValidSegCell_Trunc=ValidSegCell;
			ValidSegCell_Trunc(1,1)=num2cell(1);
			for s=2:length(OverThreshSegments)
				prevStart=ValidSegCell_Trunc(s-1,1);
				prevLength=ValidSegCell_Trunc(s-1,2);
				NewStart=prevStart{:}+prevLength{:};
				ValidSegCell_Trunc(s,1)=num2cell(NewStart);
			end	
			% save 2-column df indicating start of valid segments and length: matches dtseries at derive_wave
			segmentfnTr=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_ValidSegments_Trunc_run_' runstring];
			writetable(cell2table(ValidSegCell_Trunc),segmentfnTr,'WriteVariableNames',1)
			% make binary mask for continuous segments
			TRwise_mask_cont=zeros(1,newTRnum);
			for seg=1:length(OverThreshSegments);
				TRwise_mask_cont(ValidSegCell{seg,1}:(ValidSegCell{seg,1}+ValidSegCell{seg,2}-1))=1;
			end
			% apply to cifti
			masked_trs_cont=ts(:,logical(TRwise_mask_cont));
			newciftiSizecont=size(masked_trs_cont);
	                newTRnumcont=newciftiSizecont(2);
			% overwite diminfo of cifti
			ts_cif.diminfo{2}.length=newTRnumcont;
			% overwrite TRs for new file
			ts_cif.cdata=masked_trs_cont;
			% set output filepath
			ofp=[ResultantFolder sname '/' sname '_ses-baselineYear1Arm1_task-' task '_p2mm_masked_run_' runstring  '.dtseries.nii'];
			% There is no reason this should be a requried step
			ofp=convertStringsToChars(ofp);
			% write out motion masked cifti
			write_cifti(ts_cif,ofp);
		end
	end
end

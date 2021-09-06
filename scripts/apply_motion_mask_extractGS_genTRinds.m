function apply_motion_mask_extractGS(subj)
% this function motion masks the individual runs (non-filtered) and the concatenated (filtered)
% load in subj
topleveldir='/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/sub-*'
direc=dir(topleveldir);
% initialize empty vector for average length
TRvecNum=[];
% for each "task"
tasks=["rest","MID","SST","nback"];
for t=1:4
	task=tasks(t);
	sname=subj;
	% mask concatenated data
	fpParent=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
	fp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_bold_desc-filtered_timeseries.dtseries.nii'],'');
	% not flagging missing tasks for now, added this conditional to reflect that
	if exist(fp,'file')
		ts_cif=read_cifti(fp);
		ts=ts_cif.cdata;
		% get size
		ciftisize=size(ts);
		numTRs=ciftisize(2);
		% get cortex indices
		%CL_vertlist=ts_cif.diminfo{1}.models{1}.vertlist+1
		%CL_ind=ts_cif.diminfo{1}.models{1}.vertlist+ts_cif.diminfo{1}.models{1}.start;
                %CR_ind=ts_cif.diminfo{1}.models{2}.vertlist+ts_cif.diminfo{1}.models{2}.start;
		%C_ind=vertcat(CL_ind,CR_ind);
		C_ind=1:59412;
		% load in mask
		masfp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_desc-filteredwithoutliers_motion_mask.mat'],'');
		if exist(masfp,'file')
			mask=load(masfp);
			% get to FD_thresh of .2 mm, corresponds to threshold 21
			maskp2mm=mask.motion_data{1,21}.combined_removal;
			TRwise_mask=logical(maskp2mm);
			% length of mask corresponds to number of TRs
			% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
			TRwise_mask=~TRwise_mask;
			% remove TRs with corresp. flag
			masked_trs=ts(:,TRwise_mask);
			% reconfig cifti metadata to reflect new number of TRs
			newciftiSize=size(masked_trs);
			newTRnum=newciftiSize(2);
			% setting continuous frame threshold to 15 TRs in a row
			Threshold=15;
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
			% check that sum of TRs matches field from 3165 mask
			allRetainedSegmentTRLengths=UTSegCell(:,2);
			if (sum([allRetainedSegmentTRLengths{:}])==~mask.motion_data{1,21}.remaining_combined_count)
				error('3165 remaining combined count does not match internal representation')
			end
                        % find segments with more continuous TRs than threshold
                        OverThreshSegments=find(ContSegments>Threshold);
                        % sum remaining segments to get included TRs if this thresh chosen
                        RemainingTRs=sum(ContSegments(OverThreshSegments))
                        % index of which TR valid segments start at
                        ValidTRStarts=dInd(maskValAtChange);
			% index out segments greater than TR thresh from UnThreshSegmentCellstruct
			ValidSegCell=UTSegCell(OverThreshSegments,:);
			% adjust ValidSegCell to describe in terms of retained TRs only
			ValidSegCell(1,1)=num2cell(1);
			for s=2:length(OverThreshSegments)
				prevStart=ValidSegCell(s-1,1);
				prevLength=ValidSegCell(s-1,2);
				NewStart=prevStart{:}+prevLength{:};
				ValidSegCell(s,1)=num2cell(NewStart);
			end	
			% save 2-column df indicating start of valid segments and length
			segmentfn=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_ValidSegments'],'');
			writetable(cell2table(ValidSegCell),segmentfn,'WriteVariableNames',0)
			% overwite diminfo
			ts_cif.diminfo{2}.length=newTRnum;
			% overwrite TRs for new file
			ts_cif.cdata=masked_trs;
			% set output filepath
			ofp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_p2mm_masked.dtseries.nii'],'');
			% There is no reason this should be a requried step
			ofp=convertStringsToChars(ofp);
			% write out motion masked cifti
			write_cifti(ts_cif,ofp);
			% manually concatenate time series from individ. runs
			GSTS=zeros(1,1);
			% segment to extract and concat minimally proc. over equiv. TRs
       			WholeCortTS=zeros(59412,0);
			% for each "run", calculate global signal
			for r=1:6
				fpParent=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
				fp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_run-' string(r) '_bold_timeseries.dtseries.nii'],'');
				% not flagging missing tasks for now, added this conditional to reflect that
				if exist(fp,'file')
					ts_cif_r=read_cifti(fp);
					ts_r=ts_cif_r.cdata;
					% extract just cortex
					tsCL=ts_r(C_ind,:);
					GSTS=[GSTS mean(tsCL)];
					% segment to extract and concat minimally proc. over equiv. TRs
                                        WholeCortTS=[WholeCortTS tsCL];
				end
			end
			% remove initialization pseudovolume
			GSTS(1)=[];
			% ensure concat volumes are same size
			sizeGS=size(GSTS);
			if sizeGS(2)~=numTRs
				error('3165 and manually concatenated TR length do not match')
			end
			% use same mask as filtered/concat
			GSTS=GSTS(TRwise_mask);
			% segment to extract and concat minimally proc. over equiv. TRs
			WholeCortTS=WholeCortTS(:,TRwise_mask);
			% number of segments
        		segShape=size(ValidSegCell);
        		numSegs=segShape(1);
        		% new time series to only include valid segments
			Oords=zeros(59412,0);
			% for each segment
       			for s = 1:numSegs
        	       		% select all vertices from this segment
        	       		% first column is start of clean window, second column is length of clean window
        	     		% -1 because duration of clean window includes first TR indexed
        		        segment=WholeCortTS(:,ValidSegCell{s,1}:(ValidSegCell{s,1}+ValidSegCell{s,2}-1));
               			% append segment onto grayOrds
               			Oords = [Oords segment];
			end
			% reset cdata to be proper size - non cortical vertices set to 0
			ts_cif.cdata=zeros(91282,RemainingTRs);
			ts_cif.cdata(C_ind,:)=Oords;
        		% reset temporal dimension to reflect thresholding of short segments
        		ts_cif.diminfo{2}.length=RemainingTRs;
			% save out cifti
			WCTSfp=strjoin([fpParent,sname '_p2mm_masked_' tasks(t) '_min_proc.dtseries.nii'],'');
			write_cifti(ts_cif,char(WCTSfp));
			%%%%%%%%
			gsfp=strjoin([fpParent,sname '_p2mm_masked_' tasks(t) '_GS.csv'],'');
			writetable(array2table(GSTS),gsfp,'WriteVariableNames',0);
		else
                	missingDir=['/cbica/projects/abcdfnets/results/MissingDataReports/' sname]; 
                	mkdir(missingDir);
                	missingFile=[missingDir '/MissingData.txt'];
                	system(['echo motionMask_missing >> ' missingFile]);			
		end
	else
		missingDir=['/cbica/projects/abcdfnets/results/MissingDataReports/' sname];
		mkdir(missingDir);
		missingFile=[missingDir '/MissingData.txt'];
		system(['echo BOLD_missing >> ' missingFile]);
end
end

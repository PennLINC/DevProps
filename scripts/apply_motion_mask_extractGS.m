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
		CL_ind=ts_cif.diminfo{1}.models{1}.vertlist+ts_cif.diminfo{1}.models{1}.start;
                CR_ind=ts_cif.diminfo{1}.models{2}.vertlist+ts_cif.diminfo{1}.models{2}.start;
		C_ind=vertcat(CL_ind,CR_ind);
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
			%s-2 because we start at 3
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
			% for each "run", calculate global signal
			for r=1:5
				fpParent=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
				fp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_run-' string(r) '_bold_timeseries.dtseries.nii'],'');
				% not flagging missing tasks for now, added this conditional to reflect that
				if exist(fp,'file')
					ts_cif=read_cifti(fp);
					ts=ts_cif.cdata;
					% extract just cortex
					tsCL=ts(C_ind,:);
					GSTS=[GSTS mean(tsCL)];
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

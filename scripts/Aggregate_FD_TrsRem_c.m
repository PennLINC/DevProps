% load subjs list
Subjs=readtable('~/PWs/hcpd_subj_list.txt')

% initialize output
output=cell(height(Subjs),3);

% for each subj
for s=1:height(Subjs)
	s
	subj=table2array(Subjs(s,1));
	FDFP=['/cbica/projects/hcpd/data/motion/fd/' subj{:} '/ses-V1/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/task-carit_power_2014_FD_only.mat'];
	if isfile(FDFP)
		FD=load(FDFP);
		RemainingFD=FD.motion_data{21}.remaining_frame_mean_FD;
		output{s,1}=subj{:};
		output{s,2}=RemainingFD;
		% TRs filepath
		TRsFP=['/cbica/projects/hcpd/data/motMasked_contSegs/' subj{:} '/' subj{:} '_ses-baselineYear1Arm1_task-carit_ValidSegments_Trunc.txt'];
		if isfile(TRsFP)
			TRs=readtable(TRsFP);
			% if table is populated
			if height(TRs) > 0
				% last seg + number in last seg = number of TRs used
				lastSeg=TRs(end,1);
				TRsinlastseg=TRs(end,2);
				% -1 because second column is inclusive (starting frame counts as one in segment count)
				totTRs=table2array(lastSeg)+table2array(TRsinlastseg)-1;
				output{s,3}=totTRs;
			end
		end
	end
end
% save out
writetable(table(output),'~/PWs/Subj_FD_RemTRs_c.csv')

% load subjs list
Subjs={'sub-MSC01','sub-MSC02','sub-MSC03','sub-MSC04','sub-MSC05','sub-MSC06','sub-MSC07','sub-MSC08','sub-MSC09','sub-MSC10'};

% initialize output
output=cell(100,3);

% iterator for row placement
rowIt=0;

% for each subj
for s=1:length(Subjs)
for r=1:10
	runstrings={'01','02','03','04','05','06','07','08','09','10'};
	runstring=runstrings{r}
	s
	% make an interator so each set of values have their own row
	rowIt=rowIt+1;
	subj={table2array(Subjs(s))};
	FDFP=['/cbica/projects/pinesParcels/data/msc/' subj{:} '/ses-func' runstring '/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/task-rest_power_2014_FD_only.mat'];
	FD=load(FDFP);
	RemainingFD=FD.motion_data{41}.remaining_frame_mean_FD;
	output{rowIt,1}=[subj{:} '_run' runstring];
	output{rowIt,2}=RemainingFD;
	% TRs filepath
	TRsFP=['/cbica/projects/pinesParcels/data/msc/motMasked_contSegs/' subj{:} '/' subj{:} '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc_run_' runstring '.txt'];
	%if isfile(TRsFP)
		TRs=readtable(TRsFP);
		% if table is populated
		if height(TRs) > 0
			% last seg + number in last seg = number of TRs used
			lastSeg=TRs(end,1);
			TRsinlastseg=TRs(end,2);
			% -1 because second column is inclusive (starting frame counts as one in segment count)
			totTRs=table2array(lastSeg)+table2array(TRsinlastseg)-1;
			output{rowIt,3}=totTRs;
		end
	%end
end
end
% save out
writetable(table(output),'~/PWs/mscSubj_FD_RemTRs.csv')

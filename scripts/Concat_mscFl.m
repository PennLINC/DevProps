Subjs={'sub-MSC01','sub-MSC02','sub-MSC03','sub-MSC04','sub-MSC05','sub-MSC06','sub-MSC07','sub-MSC08','sub-MSC09','sub-MSC10'};
for s=1:length(Subjs)
	subj=Subjs{s};
	% have total length iterator for sizing output dataframe
	totesLength=0;
        % create struct for concatenated optical flow data
        us=struct();
        us.vf_left=[];
        us.vf_right=[];
	for r=1:10
		runstrings={'01','02','03','04','05','06','07','08','09','10'};
		runstring=runstrings{r}
		% load in opflow output
		OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4_run_' runstring '.mat'];
		OpFl=load(OpFlFp);
		OpFlMat=OpFl.us;
		% record duration of this run
		runDur=length(OpFlMat.vf_left);
		% yes, loop over every volumn because cells
		for c=1:runDur
			us.vf_left{c+totesLength}=OpFlMat.vf_left{c};
			us.vf_right{c+totesLength}=OpFlMat.vf_right{c};
		end
	totesLength=totesLength+runDur;
	end
	% saveout concat for PGGAngD
	OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4_concat.mat'];
	save(OpFlFp,'us')
end

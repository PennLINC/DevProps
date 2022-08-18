% just convert BUTD L and R to csv for r processing
% get subj list
subjs=readtable('/cbica/projects/pinesParcels/PWs/G300_cTRs.txt','ReadVariableNames',false);
numSubjs=height(subjs);

% for each subj
for s=1:numSubjs
        subjcell=table2array(subjs(s,1));
        subj=subjcell{:}
	% filepath 
	outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs_My_c.mat'];
	outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs_My_c.mat'];
	if isfile(outFP_L);
		% load in L
		F_L=load(outFP_L);
		% load in R
		F_R=load(outFP_R);
		% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
		matL=cell2mat(F_L.OutDf_L);
		matR=cell2mat(F_R.OutDf_R);
		csv_L=table(matL(:,[1]));
		csv_R=table(matR(:,[1]));
		writetable(csv_L,['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_My_c.csv'])
		writetable(csv_R,['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_My_c.csv'])
	end
end

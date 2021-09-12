file_pattern='/gpfs/fs001/cbica/projects/abcdfnets/results/SingleParcel_1by1/sub-NDARINV*/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/MM_vals.csv'
data_dir=dir(file_pattern);
data_locations=fullfile({data_dir.folder}, {data_dir.name});

% initialize
df=cell(length(data_locations),36);

% for each subj, extract csv into aggregate df
for i = 1:length(data_locations)
	i
	MM=readtable(char(data_locations(i)));
	MMc=table2cell(MM);
	df(i,:)=MMc;
end
writetable(cell2table(df),'~/MM_in_17Parc.csv')

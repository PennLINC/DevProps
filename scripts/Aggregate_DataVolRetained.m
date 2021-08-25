% can adapt to one of four tasks, wimpy current version needs to be executed in output Dir
filepattern='*_MID.mat' % adapt here
data_dir=dir(filepattern)   
data_locations=fullfile({data_dir.folder}, {data_dir.name});
output=cell(length(dir(filepattern)),6);   
for i=1:length(dir(filepattern)) 
data=load(char(data_locations(i)));
output(i,1:5)=data.TRcounts;
output(i,6)=data_locations(i);
end
writetable(cell2table(output),'MID-trthreshes.csv') % adapt here

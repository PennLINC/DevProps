function spin_pg(subj)
% add spin test path
addpath(genpath('/cbica/projects/abcdfnets/scripts/spin_test'));
% parent filepath to load from
sname=char(subj);
parentfp=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
% load in left and right
lPGfn=[parentfp sname '_PG_L_10k.mat'];
rPGfn=[parentfp sname '_PG_R_10k.mat'];
lPGf=load(lPGfn);
rPGf=load(rPGfn);
% convert files to just the vectors themselves
lPG=lPGf.pg_l;
rPG=rPGf.pg_r;
% write em as a table for spin test
writetable(table(lPG'),[parentfp sname 'tab_L.csv'],'WriteVariableNames',0);
writetable(table(rPG'),[parentfp sname 'tab_R.csv'],'WriteVariableNames',0);
% output of spin filename
outFn=[parentfp sname '_spunions.mat'];
% aaaand spin (x100)
SpinPermuFS([parentfp sname 'tab_L.csv'],[parentfp sname 'tab_R.csv'], 100, outFn);

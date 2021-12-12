function spin_pg(subj)
% add spin test path and cifti path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/'));
% parent filepath to load from
sname=char(subj);
parentfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];
% load in left and right
lPGfn=[parentfp sname '_PG_L_10k_rest.func.gii'];
rPGfn=[parentfp sname '_PG_R_10k_rest.func.gii'];
lPGf=gifti(lPGfn);
rPGf=gifti(rPGfn);
% convert files to just the vectors themselves
lPG=lPGf.cdata(:,1);
rPG=rPGf.cdata(:,1);
% write em as a table for spin test
writetable(table(lPG),[parentfp sname 'tab_L.csv'],'WriteVariableNames',0);
writetable(table(rPG),[parentfp sname 'tab_R.csv'],'WriteVariableNames',0);
% output of spin filename
outFn=[parentfp sname '_spunions.mat'];
% aaaand spin (x1000)
SpinPermuFS([parentfp sname 'tab_L.csv'],[parentfp sname 'tab_R.csv'], 1000, outFn);

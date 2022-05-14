addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

Subjs=readtable('~/PWs/rs_subs.csv')

for s=1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	subj=Subj{:}
	mask_mw_faces_4(subj)
end


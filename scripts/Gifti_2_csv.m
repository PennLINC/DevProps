# script name
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
L=gifti('y7_L_3k.label.gii');
R=gifti('y7_R_3k.label.gii');
writetable(table(L.cdata),'~/data/y7_L_3k.csv')
writetable(table(R.cdata),'~/data/y7_R_3k.csv')

% load in PGs, save em out as r-friendly csvs
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in euclidean distances - faces
face_ED_l=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_left_fsaverage4_faces.mat');
face_ED_r=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_right_fsaverage4_faces.mat');
% writeout
writetable(table(face_ED_l.bdsml),'~/results/PWs/face_EucD_left.csv');
writetable(table(face_ED_r.bdsmr),'~/results/PWs/face_EucD_right.csv');

% euclidean distances - vertices
eucl_l=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_left_fsaverage5.mat');
eucl_r=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_right_fsaverage5.mat');
eucl_l=eucl_l.bdsml;
eucl_r=eucl_r.bdsmr;

% writeout
writetable(table(eucl_l),'~/results/PWs/vert_EucD_left.csv');
writetable(table(eucl_r),'~/results/PWs/vert_EucD_right.csv');


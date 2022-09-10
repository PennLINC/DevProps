% in bash, used this command
% mri_annotation2label --subject fsaverage4 --hemi rh --outdir ~/data/lobes --lobes ~/data/fs4lobes_rh

% add needed paths
 addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in lobe info
[vertices_lh, label_lh, colortable_lh]=read_annotation('~/data/fs4lobes_lh.annot');
[vertices_rh, label_rh, colortable_rh]=read_annotation('~/data/fs4lobes_rh.annot');

% saveout so python can read it
dlmwrite('~/data/lobe_label_lh.csv',label_lh,'delimiter',',','precision',12);
dlmwrite('~/data/lobe_label_rh.csv',label_rh,'delimiter',',','precision',12);



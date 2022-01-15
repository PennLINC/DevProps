function PGG_AngDist(subj)

% calculate principal gradient gradient for this subject, determine angular distacne b/w PGG directionality and OpFl directionality across resting-state scan
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

%PGG_AngDistCalc(subj)

% add a script that masks medial wall w/r/t faces (not verts) out
mask_mw_faces(subj)

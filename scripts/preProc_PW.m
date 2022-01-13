function preProc_PWs(subj)

% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow
% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% apply motion masks and extract TR indcies of continuous segments
%%% already ran on subjs - NOTE HCPD DIREC VERSION RUN, EXCLUDES NONCONTINUOUS SEGMENTS 
% apply_motion_mask(subj)

%%%%%%%%%%%%%%%
% don't mess about with this stuff if it's already been ran
% fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
% PGfp=[fp '/' subj '_PG_LR_32k_rest.dscalar.nii'];
%if ~exist(PGfp)
%%%%%%%%%%%%%%%

%%% already ran on subjs
% downsample aggregated TS 
% dsCommand=['~/PWs/scripts/downsample_TS.sh ' subj];
% system(dsCommand)

% make a "Proced" dir
%direcString=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
%mkdirCommand=['mkdir ' direcString];
%system(mkdirCommand)

% derive personalized PG
% derivePGcommand=['/cbica/projects/pinesParcels/miniconda3/envs/mv_preds/bin/python derive_pg.py ' subj];
% system(derivePGcommand)

%%%%%%%%%%%%%%%
% don't mess about with OpFl if it's already been ran
OFfp=['/cbica/projects/pinesParcels/results/Proced/' subj '/' subj '_OpFl_fs5.mat'];
if ~exist(OFfp)
%%%%%%%%%%%%%%%

%%% Run Spherical Optical flow
OpFl_Sph_fs5(subj)
end

%%% Calculate the calculus gradient of the principal gradient, calculate angular distance of OpFl directions
PGG_AngDistCalc(subj)

%%% mask medial wall and extract R-friendly face data
mask_mw_faces(subj)


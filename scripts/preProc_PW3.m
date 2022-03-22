function preProc_PWs(subj)

% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow
% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% apply motion masks and extract TR indcies of continuous segments
%%% already ran on subjs - NOTE HCPD DIREC VERSION RUN, EXCLUDES NONCONTINUOUS SEGMENTS 
%apply_motion_mask(subj)

%%% already ran on subjs
% downsample aggregated TS 
%dsCommand=['~/PWs/scripts/downsample_TS.sh ' subj];
%system(dsCommand)

% make a "Proced" dir
%direcString=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
%mkdirCommand=['mkdir ' direcString];
%system(mkdirCommand)

%%%%%%%%%%%%%%%
% don't mess about with OpFl if it's already been ran
%OFfp=['/cbica/projects/pinesParcels/results/Proced/' subj '/' subj '_OpFl_fs5.mat'];
%if ~exist(OFfp)
%%%%%%%%%%%%%%%

% extract BUTD
Extract_BUTD_ResultantVecs_curv(subj)

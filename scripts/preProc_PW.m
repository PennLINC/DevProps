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

%%% RUN 1
%%% Run Spherical Optical flow
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_OpFl_Sph_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=15G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl.sh']);

% RUN 2
%%% Calculate the calculus gradient of the principal gradient, calculate angular distance of OpFl directions
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_PGG_AngDistCalc_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=13G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD.sh']);

%PGG_AngDistCalc4(subj)
%system(['/cbica/projects/pinesParcels/PWs/scripts/run_PGG_AngDistCalc_CompVer.sh $MATLAB_DIR' subj]);
%system(['/cbica/projects/pinesParcels/PWs/scripts/run_CG_AngDistCalc_CompVer.sh $MATLAB_DIR' subj]);

% RUN 3
% Spatial Nulls: PGG
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_PGG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD_snull.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=13G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD_snull.sh']);
% CG
cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_CG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR ' subj];
fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_CG_AngD_snull.sh'], 'w');
fprintf(fid,cmd);
system(['qsub -l h_vmem=15G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_CG_AngD_snull.sh']);


% RUN 4
% temporal nulls: optical flow
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_tnull_comb_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_OF_tnull.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=80G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_OF_tnull.sh']);


%%% mask medial wall and extract R-friendly face data
%mask_mw_faces_4(subj)
%mask_mw_faces_4curv(subj)

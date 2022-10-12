function preProc_PWs_msc(subj)

% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow
% print subject being ran
subj

% add matlab path for used functions
%addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% apply motion masks and extract TR indcies of continuous segments
%apply_motion_mask_msc(subj)

% downsample aggregated TS 
%dsCommand=['~/PWs/scripts/downsample_TSfs4_msc.sh ' subj];
%system(dsCommand)

% make a "Proced" dir
% direcString=['/cbica/irojects/pinesParcels/results/PWs/Proced/' subj];
% mkdirCommand=['mkdir ' direcString];
% system(mkdirCommand)


%%% RUN 1
%%% Run Spherical Optical flow
cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_OpFl_Sph_fs4_runSweep_msc.sh $MATLAB_DIR ' subj];
fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl.sh'], 'w');
fprintf(fid,cmd);
system(['qsub -l h_vmem=35G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl.sh']);

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
cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_PGG_AngDistCalc_snull_CompVer_msc.sh $MATLAB_DIR ' subj];
fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD_snull.sh'], 'w');
fprintf(fid,cmd);
system(['qsub -l h_vmem=33G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_PGG_AngD_snull.sh']);

%%%% RUN 4
%%%%% above with spins
% cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_MyG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR ' subj];
% fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=15G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG.sh']);

%%% mask medial wall and extract R-friendly face data
%mask_mw_faces_4(subj)
%mask_mw_faces_4_tertSpecific(subj)
% aggregate Bu TD prop for tert PG?

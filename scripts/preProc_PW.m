function preProc_PWs(subj)

% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow
% print subject being ran
subj

% add matlab path for used functions
%addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% apply motion masks and extract TR indcies of continuous segments
%%% already ran on subjs - NOTE HCPD DIREC VERSION RUN, EXCLUDES NONCONTINUOUS SEGMENTS 
%apply_motion_mask(subj)

%%% already ran on subjs
% downsample aggregated TS 
%dsCommand=['~/PWs/scripts/downsample_TS.sh ' subj];
%system(dsCommand)

% make a "Proced" dir
%direcString=['/cbica/irojects/pinesParcels/results/PWs/Proced/' subj];
%mkdirCommand=['mkdir ' direcString];
%system(mkdirCommand)


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
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_CG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_CG_AngD_snull.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=15G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_CG_AngD_snull.sh']);

% RUN 4
% temporal nulls: optical flow
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_tnull_comb_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_OF_tnull.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=150G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_OF_tnull.sh']);

%%% RUN 5
%%% Run Spherical Optical flow with parameter sweeps (4x the runs, intended for subset of subjs)
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_OpFl_Sph_fs4_paramSweep.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl_pS.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=25G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_OpFl_pS.sh']);

%%%% RUN 6
%%%% Run angular distance calculation on altertative reference hierarchy
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_MyG_AngDistCalc4_CompVer.sh $MATLAB_DIR ' subj];
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_MyG.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=11G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_MyG.sh']);

%%%% RUN 7
%%%%% above with spins
% cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_MyG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR ' subj];
% fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=15G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG.sh']);

%%%% RUN 8
%%%% Run angular distance calculation using participant's own tertile-derived PG
%tertiles=["young" "mid" "old"];
%tertilecodes=["Y" "Mi" "O"];
%for T=1:3
%tertile=tertiles(T)
%tertilecode=tertilecodes(T)
% load in tertile names
%subjs=readtable(['~/PWs/' tertile{:} '_subs.txt'],'ReadVariableNames',false);
%sizesubjs=size(subjs);
%lsubjs=sizesubjs(1)
% for each subject
%for s=1:lsubjs
%subj=subjs{s,1}
%cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_' tertilecode{:} 'PG_AngDistCalc4_CompVer.sh $MATLAB_DIR ' subj{:}];
% run tert ang dist calc
%fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj{:} '_tertPG.sh'], 'w');
%fprintf(fid,cmd);
%system(['qsub -l h_vmem=11G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj{:} '_tertPG.sh']);
%pause(7)
%end
%end

%%%% RUN 9
%%%% extract Myelin features

%Extract_BUTD_ResultantVecs_My(subj)

%%% RUN 10
% run angular distance on myelin over tasks
cmd=['/cbica/projects/pinesParcels/PWs/scripts/run_MyG_AngDistCalc_c_CompVer.sh $MATLAB_DIR ' subj];
fid=fopen(['/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG_c.sh'], 'w');
fprintf(fid,cmd);
system(['qsub -l h_vmem=13G ' '/cbica/projects/pinesParcels/data/CombinedData/' subj '_SMyG_c.sh']);


%%% run 11:
%%% mask medial wall and extract R-friendly face data
%mask_mw_faces_4(subj)

% aggregate Bu TD prop for tert PG?

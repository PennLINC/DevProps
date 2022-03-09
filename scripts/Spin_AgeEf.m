% addpaths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
% set outdir
outdir='/cbica/projects/pinesParcels/results/aggregated_data/';

% set up masking indices
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[V_L, F_L] = read_surf(surfL);
[V_R, F_R] = read_surf(surfR);
% +1 the faces: begins indexing at 0
F_L = F_L + 1;
F_R = F_R + 1;
% load in PG, PGG it for masking
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

% get Age effect map
ageEfL = readtable(['/cbica/projects/pinesParcels/results/PWs/FDRed_Prop_L.csv'],'ReadVariableNames',0);
ageEfR = readtable(['/cbica/projects/pinesParcels/results/PWs/FDRed_Prop_R.csv'],'ReadVariableNames',0);
% plop into non-mw corresponding spots
dataL=zeros(1,5120);
dataL(g_noMW_combined_L)=table2array(ageEfL);
dataR=zeros(1,5120);
dataR(g_noMW_combined_R)=table2array(ageEfR);

% initialize vertVecs
vertVecL=zeros(1,length(V_L));
vertVecR=zeros(1,length(V_R));

% convert faces to vertices
for vL=1:length(V_L)
	[rowInd,cInd]=find(F_L==vL);
	vertVecL(vL)=mean(dataL(rowInd));	
end
% right hemi
for vR=1:length(V_R)
	[rowInd,cInd]=find(F_R==vR);
        vertVecR(vR)=mean(dataR(rowInd));
end

% set masked regions to 1000 (using vertexwise mw mask)
surfML = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
vertVecL(mwIndVec_l)=1000;
vertVecR(mwIndVec_r)=1000;

% set output file name
outFn=strcat('/cbica/projects/pinesParcels/results/aggregated_data/AgeEf_Permuts_fs4.mat');

% write them out as a csv for spin test to deal with		
writetable(table(vertVecL),[outdir 'AgeEfs_fs4_L.csv'],'WriteVariableNames',0);
writetable(table(vertVecR),[outdir 'AgeEfs_fs4_R.csv'],'WriteVariableNames',0);
% create permutations, save out to outFn
SpinPermuFS([outdir 'AgeEfs_fs4_L.csv'], [outdir 'AgeEfs_fs4_R.csv'], 1000, outFn,'_AgeEf_');

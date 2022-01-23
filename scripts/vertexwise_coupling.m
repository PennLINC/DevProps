%function facewise_coupling(subj)

addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

%%% load in subject ts
% load in TRs_l and TRs_r
TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_10k.mgh'];
TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_10k.mgh'];
% filepaths to files
TRs_lf=MRIread(TRs_lfp);
TRs_rf=MRIread(TRs_rfp);
% files to data
TRs_l=squeeze(TRs_lf.vol);
TRs_r=squeeze(TRs_rf.vol);

%%%% load in fsaverage5 faces and vertices
% for surface data
SubjectsFolder='/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));

%%% load in Angular distance from PGG
angDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat.mat'];
Angs=load(angDistFP);
AngsL=Angs.AngDist.gLeft';
AngsR=Angs.AngDist.gRight';

%%% NEED TO CONVERT ANG DIST TO VERTS RATHER THAN TS TO FACES

% time series size
sizeTsl=size(TRs_l);
lengthTSl=sizeTsl(2);
sizeTsr=size(TRs_r);
lengthTSr=sizeTsr(2);
% and for opflow
OsizeTsl=size(AngsL);
OlengthTSl=OsizeTsl(2);
OsizeTsr=size(AngsR);
OlengthTSr=OsizeTsr(2);

% initialize angular vectors
VertAngsL=zeros(length(vx_l),lengthTSl);
VertAngsR=zeros(length(vx_r),lengthTSr);

for lv=1:length(vx_l)
	% get instances where this vert is part of face column 1
	inds1=find(faces_l(:,1)==lv);
        inds2=find(faces_l(:,2)==lv);
        inds3=find(faces_l(:,3)==lv);
	allinds=union(inds1,inds2);
	allinds=union(allinds,inds3);
	% average value of involved faces at each TR
	for TR=1:OlengthTSl
		VertAngsL(lv,TR)=mean(AngsL(allinds,TR));
	end
end

for rv=1:length(vx_r)
        % get instances where this vert is part of face column 1
        inds1=find(faces_r(:,1)==rv);
        inds2=find(faces_r(:,2)==rv);
        inds3=find(faces_r(:,3)==rv);
        allinds=union(inds1,inds2);
        allinds=union(allinds,inds3);
	for TR=1:OlengthTSr
        	VertAngsR(rv,TR)=mean(AngsR(allinds,TR));
	end
end

% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);
% MASK WHERE PGG = 0: individ AND group
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.lh.fsaverage5.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.rh.fsaverage5.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gLPGf.cdata(:,1:10)')==0);
gPGg_R0=find(all(gRPGf.cdata(:,1:10)')==0);
% get unions of MW masks and empty faces
gro_MW_combined_L=union(mw_v_l,gPGg_L0);
gro_MW_combined_R=union(mw_v_r,gPGg_R0);
% pull in individualized PGs to make sure we are masking comprehensively
% load in subject's PG
LPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_L_10k_rest.func.gii'];
LPGf=gifti(LPGfp);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_R_10k_rest.func.gii'];
RPGf=gifti(RPGfp);
PGg_L0=find(all(LPGf.cdata(:,1:10)')==0);
PGg_R0=find(all(RPGf.cdata(:,1:10)')==0);
% combined gorup and individ and mw masks
MW_combined_L=union(gro_MW_combined_L,PGg_L0);
MW_combined_R=union(gro_MW_combined_R,PGg_R0);
% get opposite: where we don't want to mask
noMW_combined_L=setdiff([1:10242],MW_combined_L);
noMW_combined_R=setdiff([1:10242],MW_combined_R);



% Mask vert TS
vertTS_l=TRs_l(noMW_combined_L,:);
vertTS_r=TRs_r(noMW_combined_R,:);
% mask hierarchical angling
AngsL=VertAngsL(noMW_combined_L,:);
AngsR=VertAngsR(noMW_combined_R,:);
display('Done masking, now computing')

%%%% compute FC matrix
bothHemisTS=vertcat(vertTS_l,vertTS_r);
VertFC=corrcoef(bothHemisTS');

% compute hierarchical anglular coupling
bothHemisHierAngTS=vertcat(AngsL,AngsR);
AngFC=corrcoef(bothHemisHierAngTS');
% make an upper triangle mask to avoid diagonals and redundancies
TriuMask=triu(ones(18594,18594),1);
% vectors of both
VertFCvec=VertFC(logical(TriuMask));
AngFCvec=AngFC(logical(TriuMask));







%%%%%%%%%%%%%%%%%%%% Distance vectors

gPHmerged=vertcat(gPG_LH(noMW_combined_L),gPG_RH(noMW_combined_R));
pgDistmat=zeros(18594,18594);
for i=1:18594
	for j=1:18594
		pgDistmat(i,j)=abs(gPHmerged(i)-gPHmerged(j));
	end
end
% same upper triangle mask
pgdistvec=pgDistmat(logical(TriuMask));


%%%% Convert Angular FC and FC back to single-hemisphere for comparison
LHfc=corrcoef(vertTS_l');
RHfc=corrcoef(vertTS_r');


% now get a distance vector
eucl_l=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_left_fsaverage5.mat');
eucl_r=load('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_right_fsaverage5.mat');
eucl_l=eucl_l.bdsml;
eucl_r=eucl_r.bdsmr;

% import connectome-workbench-generated geodesic distance matrices (calc on fsaverage5.surf.gii)
lhfile=ciftiopen('/cbica/projects/pinesParcels/data/lh_GeoDist.dconn.nii','/cbica/software/external/connectome_workbench/1.4.2/bin/wb_command');
rhfile=ciftiopen('/cbica/projects/pinesParcels/data/rh_GeoDist.dconn.nii','/cbica/software/external/connectome_workbench/1.4.2/bin/wb_command');
lhGeoDmat=lhfile.cdata;
rhGeoDmat=rhfile.cdata;


CouplingTable=table(VertFCvec,AngFCvec);

% save out both as tall csv for this subj (for R)
writetable(CouplingTable,['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_CouplingTable.csv']);

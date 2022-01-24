function facewise_coupling(subj)

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
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5';
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

%%%% convert vertex ts to face ts
sizeTsl=size(TRs_l);
lengthTSl=sizeTsl(2);
sizeTsr=size(TRs_r);
lengthTSr=sizeTsr(2);
% initialize face TS : # of faces and # of TRs
faceTS_l=zeros(length(faces_l),sizeTsl(2));
faceTS_r=zeros(length(faces_r),sizeTsr(2));
% there is probably a more efficient way to do this than loop over every TR
for t=1:lengthTSl
	ts_this_tr_l=TRs_l(:,t);
	faceTS_l(:,t)=sum(ts_this_tr_l(faces_l),2)./3;
        ts_this_tr_r=TRs_r(:,t);
        faceTS_r(:,t)=sum(ts_this_tr_r(faces_r),2)./3;
end

%%% make like Jacques Plante and bring in face masks
% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);
% all faces that touch a medial wall vertex to be masked
MW_f1_L=find(ismember(faces_l(:,1),mw_v_l));
MW_f2_L=find(ismember(faces_l(:,2),mw_v_l));
MW_f3_L=find(ismember(faces_l(:,3),mw_v_l));
% rh
MW_f1_R=find(ismember(faces_r(:,1),mw_v_r));
MW_f2_R=find(ismember(faces_r(:,2),mw_v_r));
MW_f3_R=find(ismember(faces_r(:,3),mw_v_r));
% inclusive - mask if face involves ANY mw vertices
MW_combined_L=union(MW_f1_L,MW_f2_L);
MW_combined_L=union(MW_combined_L,MW_f3_L);
% now for right hemisphere
MW_combined_R=union(MW_f1_R,MW_f2_R);
MW_combined_R=union(MW_combined_R,MW_f3_R);
% MASK WHERE PGG = 0: individ AND group
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.lh.fsaverage5.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.rh.fsaverage5.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% calculate group PG gradient on sphere
gPGg_L = grad(faces_l, vx_l, gPG_LH);
gPGg_R = grad(faces_r, vx_r, gPG_RH);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
% get unions of MW masks and empty faces
gro_MW_combined_L=union(MW_combined_L,gPGg_L0);
gro_MW_combined_R=union(MW_combined_R,gPGg_R0);
% pull in individualized PGs to make sure we are masking comprehensively
% load in subject's PG
LPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_L_10k_rest.func.gii'];
LPGf=gifti(LPGfp);
PG_LH=LPGf.cdata(:,1);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_R_10k_rest.func.gii'];
RPGf=gifti(RPGfp);
PG_RH=RPGf.cdata(:,1);
% calculate PG gradient on sphere
PGg_L = grad(faces_l, vx_l, PG_LH);
PGg_R = grad(faces_r, vx_r, PG_RH);
PGg_L0=find(all(PGg_L')==0);
PGg_R0=find(all(PGg_R')==0);
% combined gorup and individ and mw masks
MW_combined_L=union(gro_MW_combined_L,PGg_L0);
MW_combined_R=union(gro_MW_combined_R,PGg_R0);
% get opposite: where we don't want to mask
noMW_combined_L=setdiff([1:20480],MW_combined_L);
noMW_combined_R=setdiff([1:20480],MW_combined_R);
% Mask face TS
faceTS_l=faceTS_l(noMW_combined_L,:);
faceTS_r=faceTS_r(noMW_combined_R,:);
% mask hierarchical angling
AngsL=AngsL(noMW_combined_L,:);
AngsR=AngsR(noMW_combined_R,:);
display('Done masking, now computing')

%%%% compute FC matrix
bothHemisTS=vertcat(faceTS_l,faceTS_r);
FaceFC=corrcoef(bothHemisTS');
% compute hierarchical anglular coupling
bothHemisHierAngTS=vertcat(AngsL,AngsR);
AngFC=corrcoef(bothHemisHierAngTS');
% make an upper triangle mask to avoid diagonals and redundancies
TriuMask=triu(ones(37059,37059),1);
% vectors of both
FaceFCvec=FaceFC(logical(TriuMask));
AngFCvec=AngFC(logical(TriuMask));
CouplingTable=table(FaceFCvec,AngFCvec);

% save out both as tall csv for this subj (for R)
writetable(CouplingTable,['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_CouplingTable.csv']);

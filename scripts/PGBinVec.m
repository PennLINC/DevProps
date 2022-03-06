function PGBinVec(VecL,VecR,fnout)
% take vectors and save out csvs for each pgbin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load in gPPG ang dist: for thresholdingdd OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% freesurfer info
% Load in surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);

%%% For masking out the medial wall
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
% get inverse for indexing : faces that ARE NOT touching mW verts
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

% convert vertices to faces
PGL=sum(gPG_LH(faces_l),2)./3;
PGR=sum(gPG_RH(faces_r),2)./3;

% combine PG for perc binning
PG=vertcat(PGL(g_noMW_combined_L),PGR(g_noMW_combined_R));
Vec=horzcat(VecL,VecR);
% for each 10% stretch
prctiles=prctile(PG,(0:20:100));

% for each percentile
for p=1:10
	% get percentile bounds
	boundl=prctiles(p);
	boundu=prctiles(p+1);
	% get index of where in PG these bounds are met
	idx = find(PG<boundu & PG > boundl);
	% extract bin values from VecL and VecR at those indices
	BinVals=Vec(idx);
	% saveout as fnout _p.csv
	fn=[fnout num2str(p) '.csv'];
	writetable(table(BinVals),fn);
end

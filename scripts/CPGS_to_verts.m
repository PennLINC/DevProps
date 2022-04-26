% load left joel embeds
embL=load('/cbica/projects/pinesParcels/results/PWs/embL.mat');
% right
embR=load('/cbica/projects/pinesParcels/results/PWs/embR.mat');

% % Load in surface data
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

%%% embeddings are in non-medial wall indices, load in medial-wall to back-construct
% use native freesurfer command for mw mask indices
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% initialize "full" vectors for first 5 PGs
vecsL=zeros(5,5120);
vecsR=zeros(5,5120);
vecL(1:5,g_noMW_combined_L)=embL.emb(:,1:5)';
vecR(1:5,g_noMW_combined_R)=embR.emb(:,1:5)';

% initialize vertexwise vectors
LoutputVerts=zeros(5,2562);
RoutputVerts=zeros(5,2562);

% for each PG
for p=1:5
	for v=1:2562
		[rowL,colL]=find(faces_l==v);
		vertVals=mean(vecL(rowL));
		LoutputVerts(p,v)=vertVals;
	end
	for v=1:2562
                [rowR,colR]=find(faces_r==v);
                vertVals=mean(vecR(rowR));
                RoutputVerts(p,v)=vertVals;
        end
end

writetable(table(LoutputVerts),'Left_circPGs.csv');
writetable(table(RoutputVerts),'Right_circPGs.csv');




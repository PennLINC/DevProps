Modes=readtable('test_leftmode.csv')


addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in angular distances
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat_L.mat'];
data=load(AngDistFP)

%%% Load in surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5';
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

% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);

% all faces that touch a medial wall vertex to be masked
MW_f1_L=find(ismember(F_L(:,1),mw_v_l));
MW_f2_L=find(ismember(F_L(:,2),mw_v_l));
MW_f3_L=find(ismember(F_L(:,3),mw_v_l));
% rh
MW_f1_R=find(ismember(F_R(:,1),mw_v_r));
MW_f2_R=find(ismember(F_R(:,2),mw_v_r));
MW_f3_R=find(ismember(F_R(:,3),mw_v_r));
% inclusive - mask if face involves ANY mw vertices
MW_combined_L=union(MW_f1_L,MW_f2_L);
MW_combined_L=union(MW_combined_L,MW_f3_L);
% now for right hemisphere
MW_combined_R=union(MW_f1_R,MW_f2_R);
MW_combined_R=union(MW_combined_R,MW_f3_R);
% get inverse for indexing : faces that ARE NOT touching mW verts
noMW_combined_L=setdiff([1:20480],MW_combined_L);
noMW_combined_R=setdiff([1:20480],MW_combined_R);









mycolormap = customcolormap([0 .2 .4 .6 .8 1], {'#ff0000','#ffa500','#ffff00','#00ff00','#00ffff','#0000ff'});


%%%%
data=zeros(1,20480);
data(noMW_combined_L)=table2array(Modes(:,2));



[vertices, faces] = freesurfer_read_surf('/cbica/software/external/freesurfer/scientificlinux6/6.0.0/subjects/fsaverage5/surf/lh.inflated');


figure
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([270 0]);
colormap(mycolormap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong;
shading flat;
camlight;
alpha(1)
colorbar




set(gca,'CLim',[0,24]);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


print('test_subj2.png','-dpng')

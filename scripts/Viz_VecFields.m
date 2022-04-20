%%%%%% Set paths %%%
% grab tool dir
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
% set fs dir
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
% set output dir
vizdir='/cbica/projects/pinesParcels/results/viz/waves/';
% subtight plot - thank you Felipe G Nievinski!
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end
% https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set "subj" variable manually: subj='sub-xxx'
%%% load op flo from subject %%%
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
OpFl=load(OpFlFp);
OpFl=OpFl.us;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Load TS %%%%%%%%
% load in TRs_l and TRs_r
TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_3k.mgh'];
TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_3k.mgh'];
% filepaths to files
TRs_lf=MRIread(TRs_lfp);
TRs_rf=MRIread(TRs_rfp);
% files to data
TRs_l=squeeze(TRs_lf.vol);
TRs_r=squeeze(TRs_rf.vol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% uptake functional data (on surface) %%%%%%
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l)
sizeInDr=size(TRs_r)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')
% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l(:,TRP);
end
% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r(:,TRP);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in fsaverage4 faces and vertices %%%
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% Get incenters of triangles.
TR = TriRep(faces_l, vx_l);
P = TR.incenters;
TRr = TriRep(faces_r, vx_r);
Pr = TRr.incenters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load in group PG %%%%%
LPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
LPGf=gifti(LPGfp);
PG_LH=LPGf.cdata(:,1);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
RPGf=gifti(RPGfp);
PG_RH=RPGf.cdata(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% mask angular distances
goodIndr=g_noMW_combined_R;
% calculate PG gradient on sphere
PGg_L = grad(faces_l, vx_l, PG_LH);
PGg_R = grad(faces_r, vx_r, PG_RH);
% extract face-wise vector cartesian vector components
PGx_L=PGg_L(:,1);
PGy_L=PGg_L(:,2);
PGz_L=PGg_L(:,3);
PGx_R=PGg_R(:,1);
PGy_R=PGg_R(:,2);
PGz_R=PGg_R(:,3);
%%% set colormap
custommap=colormap('inferno'); %or whatever
% reduce just a little bit on the close-to-white coloring
custommap=custommap(1:240,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Fig 1a: example optical flow vectors
% roy-big-bl palette imitation, inferno is just template
roybigbl_cm=inferno(16);
roybigbl_cm(1,:)=[255, 255, 0 ];
roybigbl_cm(2,:)=[255, 200, 0];
roybigbl_cm(3,:)=[255, 120, 0];
roybigbl_cm(4,:)=[255, 0, 0 ];
roybigbl_cm(5,:)=[200, 0, 0 ];
roybigbl_cm(6,:)=[150, 0, 0 ];
roybigbl_cm(7,:)=[100, 0, 0 ];
roybigbl_cm(8,:)=[60, 0, 0 ];
roybigbl_cm(9,:)=[0, 0, 80 ];
roybigbl_cm(10,:)=[0, 0, 170];
roybigbl_cm(11,:)=[75, 0, 125];
roybigbl_cm(12,:)=[125, 0, 160];
roybigbl_cm(13,:)=[75, 125, 0];
roybigbl_cm(14,:)=[0, 200, 0];
roybigbl_cm(15,:)=[0, 255, 0];
roybigbl_cm(16,:)=[0, 255, 255];
% scale to out of 1
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate in-between values
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
roybigbl_cm=roybigbl_cm(15:240,:);
% pulled from https://github.com/Washington-University/workbench/blob/master/src/Files/PaletteFile.cxx

%%% Load in subject time series indices
parentfp=['/cbica/projects/hcpd/data/motMasked_contSegs/'];
CSIfp = [parentfp subj '/' subj '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'];
CSI=readtable(CSIfp);
% get index of last TR in cont seg
lastInSegs=CSI{:,1}+CSI{:,2}-1;
% set diff between these lastInSegs and sequence of 1:#trs (-1 because of inclusivity of second column)
% go to motion masking scripts for more detail on that
numTrs=CSI{end,1}+CSI{end,2}-1;
% invalid TR pairs are those after the last TR in segments
validTRs=setdiff([1:numTrs],lastInSegs);
% now we should be able to index the desired TR based on the tr pair
for i=444:452
OpFlVecofInt=i;
TRofInt=validTRs(OpFlVecofInt);
u=OpFl.vf_right{OpFlVecofInt};
vATTR=fr.TRs{TRofInt};
% z-score
vATTR=zscore(vATTR);
% normalize vectors, pull em in
figure('units','pixels','position',[0 0 4500 4500])
axis([-1, 1, -1, 1, 0, 1]);
% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=5;
% convert vectors to unit vectors for plotting independent of magnitude
% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
ret = bsxfun(@rdivide, u, sqrt(sum(u'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3)],'w',.7,'arrowRadius',.05)
hold on
% for OpFl Vecs on PG
%trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), PG_RH, 'EdgeColor','none');
%caxis([-5.5,6.5]);
% for OpFl Vecs on BOLD
trisurf(faces_r, vx_r(:, 1)./scalingfactor, vx_r(:, 2)./scalingfactor, vx_r(:, 3)./scalingfactor, vATTR, 'EdgeColor','none');
caxis([-3,3])
axis equal
daspect([1, 1, 1]);
colormap(roybigbl_cm);
c=colorbar
c.FontSize=55
c.LineWidth=3
c.Ticks=[-3 -2 -1 0 1 2 3]
c.Location='southoutside'
c.FontName='Arial'
view(280,185);
fn=['yourfigure' num2str(i) '.png'];
print(fn,'-dpng')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Fig 1b: PGG big
custommap=colormap('inferno'); %or whatever
% reduce just a little bit on the close-to-white coloring
custommap=custommap(1:240,:);
figure('units','pixels','position',[0 0 4500 4500])
axis([-1, 1, -1, 1, 0, 1]);
% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=5;
% convert vectors to unit vectors for plotting independent of magnitude
% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
ret = bsxfun(@rdivide, PGg_R, sqrt(sum(PGg_R'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3)],'w',.7,'arrowRadius',.05)
%quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), PGx_R, PGy_R, PGz_R, 2, 'w','linewidth',2);
hold on
trisurf(faces_r, vx_r(:, 1)./scalingfactor, vx_r(:, 2)./scalingfactor, vx_r(:, 3)./scalingfactor, PG_RH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colormap(custommap);
colorbar
view(280,185);
% might take > 10 minutes to render
print('pggBig.png','-dpng')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2b: BOLD vectors on PGG for Bottom up
custommap=colormap('inferno'); %or whatever
% reduce just a little bit on the close-to-white coloring
custommap=custommap(1:240,:);
parentfp=['/cbica/projects/hcpd/data/motMasked_contSegs/'];
CSIfp = [parentfp subj '/' subj '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'];
CSI=readtable(CSIfp);
% get index of last TR in cont seg
lastInSegs=CSI{:,1}+CSI{:,2}-1;
% set diff between these lastInSegs and sequence of 1:#trs (-1 because of inclusivity of second column)
% go to motion masking scripts for more detail on that
numTrs=CSI{end,1}+CSI{end,2}-1;
% invalid TR pairs are those after the last TR in segments
validTRs=setdiff([1:numTrs],lastInSegs);
for i=217:219
OpFlVecofInt=i;
TRofInt=validTRs(OpFlVecofInt);
u=OpFl.vf_right{OpFlVecofInt};
vATTR=fr.TRs{TRofInt};
% z-score
vATTR=zscore(vATTR);
% normalize vectors, pull em in
figure('units','pixels','position',[0 0 4500 4500])
axis([-1, 1, -1, 1, 0, 1]);
% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=5;
% convert vectors to unit vectors for plotting independent of magnitude
% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
ret = bsxfun(@rdivide, u, sqrt(sum(u'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3)],'w',.7,'arrowRadius',.05)
hold on
% for OpFl Vecs on PG
trisurf(faces_r, vx_r(:, 1)./scalingfactor, vx_r(:, 2)./scalingfactor, vx_r(:, 3)./scalingfactor, PG_RH, 'EdgeColor','none');
caxis([-5.5,6.5]);
axis equal
daspect([1, 1, 1]);
colormap(custommap);
c=colorbar
c.FontSize=55
c.LineWidth=3
c.Location='southoutside'
c.FontName='Arial'
view(280,185);
fn=['yourfigure' num2str(i) '.png'];
print(fn,'-dpng')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2c: BOLD vectors on PGG for Top Down
custommap=colormap('inferno'); %or whatever
% reduce just a little bit on the close-to-white coloring
custommap=custommap(1:240,:);
parentfp=['/cbica/projects/hcpd/data/motMasked_contSegs/'];
CSIfp = [parentfp subj '/' subj '_ses-baselineYear1Arm1_task-rest_ValidSegments_Trunc.txt'];
CSI=readtable(CSIfp);
% get index of last TR in cont seg
lastInSegs=CSI{:,1}+CSI{:,2}-1;
% set diff between these lastInSegs and sequence of 1:#trs (-1 because of inclusivity of second column)
% go to motion masking scripts for more detail on that
numTrs=CSI{end,1}+CSI{end,2}-1;
% invalid TR pairs are those after the last TR in segments
validTRs=setdiff([1:numTrs],lastInSegs);
for i=1058:1060
OpFlVecofInt=i;
TRofInt=validTRs(OpFlVecofInt);
u=OpFl.vf_right{OpFlVecofInt};
vATTR=fr.TRs{TRofInt};
% z-score
vATTR=zscore(vATTR);
% normalize vectors, pull em in
figure('units','pixels','position',[0 0 4500 4500])
axis([-1, 1, -1, 1, 0, 1]);
% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=5;
% convert vectors to unit vectors for plotting independent of magnitude
% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
ret = bsxfun(@rdivide, u, sqrt(sum(u'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3)],'w',.7,'arrowRadius',.05)
hold on
% for OpFl Vecs on PG
trisurf(faces_r, vx_r(:, 1)./scalingfactor, vx_r(:, 2)./scalingfactor, vx_r(:, 3)./scalingfactor, PG_RH, 'EdgeColor','none');
caxis([-5.5,6.5]);
axis equal
daspect([1, 1, 1]);
colormap(custommap);
c=colorbar;
c.FontSize=55
c.LineWidth=3
c.Location='southoutside'
c.FontName='Arial'
view(280,185);
fn=['yourfigure' num2str(i) '.png'];
print(fn,'-dpng')
end



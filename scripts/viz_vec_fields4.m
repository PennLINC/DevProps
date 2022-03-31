function viz_vec_fields4(subj)
% visualize calculated vector fields on fsaverage4 sphere

% 1 is vector fields overlaid onto spherical harmonic time series representation, left and righ hemisphere

% 2 is PGG, mag and dir

% 3 is a grayOrd plot: vertices (arranged by PG) over time

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
% and individual
iLPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_L_10k_rest.func.gii'];
iLPGf=gifti(iLPGfp);
iPG_LH=iLPGf.cdata(:,1);
% right hemi
iRPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_R_10k_rest.func.gii'];
iRPGf=gifti(iRPGfp);
iPG_RH=iRPGf.cdata(:,1);

% get sorting order of PG
[outL,idxL] = sort(iPG_LH);
[outR, idxR] = sort(iPG_RH);

% mw index
mwIndr=find(outR==0);
goodIndr=setdiff([1:length(outR)],mwIndr);

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

%%%% Fig 1: PGG big
figure('units','pixels','position',[0 0 1500 1500])
axis([-1, 1, -1, 1, 0, 1]);
quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), PGx_R, PGy_R, PGz_R, 2, 'w','linewidth',2);
hold on
trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), PG_RH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colormap(custommap);
colorbar
view(280,185);
print('pggBig.png','-dpng')

%%% Fig 1: pial PG
% load in fsaverage4 faces and vertices %%%
surfL = [SubjectsFolder '/surf/lh.pial'];
surfR = [SubjectsFolder '/surf/rh.pial'];
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
figure('units','pixels','position',[0 0 1500 1500])
axis([-1, 1, -1, 1, 0, 1]);
hold on
trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), PG_RH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colorbar
colormap(custommap)
view(280,185);
print('pggPial.png','-dpng')

% Fig S2: CG Big
cg=load('~/data/fs4curv.mat');
cG_LH=cg.gpg.gPG_LH;
cG_RH=cg.gpg.gPG_RH;
% calculate group curvature gradient on sphere
cg_L = grad(faces_l, vx_l, cG_LH);
cg_R = grad(faces_r, vx_r, cG_RH);
% extract face-wise vector cartesian vector components
cGx_L=cg_L(:,1);
cGy_L=cg_L(:,2);
cGz_L=cg_L(:,3);
cGx_R=cg_R(:,1);
cGy_R=cg_R(:,2);
cGz_R=cg_R(:,3);

figure('units','pixels','position',[0 0 1500 1500])
axis([-1, 1, -1, 1, 0, 1]);
quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), cGx_R, cGy_R, cGz_R, 2, 'w','linewidth',2);
hold on
trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), cG_RH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colormap(custommap);
colorbar
view(280,185);
print('cgBig.png','-dpng')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%%%% This chunk is for video: not needed for figs
% vector field
DirecsVecs=struct('cdata',[],'colormap',[]);
speshframes=[209:226];
for i=1:17
j=speshframes(i)
u=OpFl.vf_right{j};
vATTR=fr.TRs{j};
vATTR=zscore(vATTR);
%%%%%%
figure('units','pixels','position',[0 0 1500 1500])
subplot(2,2,1)
axis([-1, 1, -1, 1, 0, 1]);
quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), u(:, 1), u(:, 2), u(:, 3), 2, 'w','linewidth','3');
hold on
trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), vATTR, 'EdgeColor','none');
axis equalc
daspect([1, 1, 1]);
caxis([-3,3]);
colorbar
view(270,200);
%%%colormap(custommap)
%%%%
colorbar
% medial
%view(60,190)
% this is going to be easier to just photoshop out the axes if needed - appears to remove quiver3 coloring
%axis('off')
% lateral?
view(280,185);
% 3
subplot(2,2,[2 4])
limz=[-100 100];
% get time series org. by PG
PGts=TRs_r(idxR,:);
% extract segment of PGts
PGts=PGts(goodIndr,180:230);
imagesc(PGts,limz);
hold on;
% add formal lookup with pair-TR correspondence established below
line([224,224], [0,10000], 'Color', 'w');
DirecsVecs(i)=getframe(gcf);
end

% create videowriter object
video = VideoWriter([vizdir 'testDirecVecs_fs4.avi'],'Uncompressed AVI');
video.FrameRate = 1;

% open it, plop Direcs in
open(video)
writeVideo(video, DirecsVecs);
close(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% pulled from https://github.com/Washington-University/workbench/blob/master/src/Files/PaletteFile.cxx
% scale to 1
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
roybigbl_cm=roybigbl_cm(15:240,:);
colormap(roybigbl_cm);
print('yourfigure.png','-dpng')


%%%%% create conversion vector: the aim is to go from included TRs to included TR pairs (valid adjacent frames).
% the main discrep. comes from discontinuous segments: if frames 1:5 are really 1,2,4,5,6, (3 had motion)
% then  pairs (1,2) (4,5) and (5,6) are analyzed. 
% so pull out non-valid TR pairs (i.e., (2,3)) and setdiff to get index of OpFl estimations w/r/t retained TRs
% the last TR of a continuous segments does not have an opfl vec field ascribed to it
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
figure('units','pixels','position',[0 0 800 800])
axis([-1, 1, -1, 1, 0, 1]);
quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), u(:, 1), u(:, 2), u(:, 3), 2, 'w');
hold on
% for OpFl Vecs on PG
%trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), PG_RH, 'EdgeColor','none');
%caxis([-5.5,6.5]);
% for OpFl Vecs on BOLD
trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), vATTR, 'EdgeColor','none');
caxis([-3,3])
axis equal
daspect([1, 1, 1]);
%colormap(roybigbl_cm);
colormap(custommap);
c=colorbar
c.FontSize=55
c.LineWidth=3
c.Ticks=[-3 -2 -1 0 1 2 3]
c.Location='southoutside'
c.FontName='Arial'
view(280,185);
%view(60,190)
fn=['yourfigure' num2str(i) '.png'];
print(fn,'-dpng')
end

%%% and carpet plots redo
figure('units','pixels','position',[0 0 5600 2600])
limz=[-100 100];
% get time series org. by PG
PGts=TRs_r(idxR,:);
%PGts=zscore(PGts);
imagesc(PGts,limz);
colormap(roybigbl_cm);
%hold on;
%line([j,j+2], [0,10000], 'Color', 'r');
print('carpet2.png','-dpng')

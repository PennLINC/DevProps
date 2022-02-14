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
% load in fsaverage5 faces and vertices %%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load in group PG %%%%%
LPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
LPGf=gifti(LPGfp);
PG_LH=LPGf.cdata(:,1);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
RPGf=gifti(RPGfp);
PG_RH=RPGf.cdata(:,1);

% get sorting order of PG
[outL,idxL] = sort(PG_LH);

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

% 1 - VECTOR FIELDS OVERLAID ONTO ORIGINAL TRS
% vector field
DirecsVecs=struct('cdata',[],'colormap',[]);
for j=1:100
j
u=OpFl.vf_left{j};
vATTR=fl.TRs{j};
figure('units','pixels','position',[0 0 1000 1000])
subplot(2,2,1)
axis([-1, 1, -1, 1, 0, 1]);
quiver3(P(:, 1), P(:, 2), P(:, 3), u(:, 1), u(:, 2), u(:, 3), 3, 'k');
hold on
trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), vATTR, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
caxis([-45,45]);
colorbar
view(200,200);

% 2
subplot(2,2,3)
axis([-1, 1, -1, 1, 0, 1]);
quiver3(P(:, 1), P(:, 2), P(:, 3), PGx_L, PGy_L, PGz_L, 2, 'w');
hold on
trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), PG_LH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colormap(custommap)
colorbar
view(225,215);

% 3
subplot(2,2,[2 4])
limz=[-100 100];
% get time series org. by PG
PGts=TRs_l(idxL,:);
imagesc(PGts,limz);
hold on;
line([j,j+2], [0,10000], 'Color', 'r');


DirecsVecs(j)=getframe(gcf);
end

% create videowriter object
video = VideoWriter([vizdir 'testDirecVecs_fs4.avi'],'Uncompressed AVI');
video.FrameRate = 2;

% open it, plop Direcs in
open(video)
writeVideo(video, DirecsVecs);
close(video);


function viz_vec_fields4(subj)
% visualize calculated vector fields on fsaverage4 sphere

% 1 is vector fields overlaid onto spherical harmonic time series representation, left and righ hemisphere

% 2 is PGG, mag and dir


%%%%%% Set paths %%%
% grab tool dir
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
% set fs dir
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';

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
%%%%%%% uptake functi nal data (on surface) %%%%%%
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

%%%% faces-to-vertices converter

% init vertices cell array with []s for each vertex as a cell
% for each face
% extract vertices comrpising F
% append vertex1 cell 
% append vertex2 cell
% append vertex3 cell
% end
% average vector within each vertex cell for vertex value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load in group PG %%%%%
LPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
LPGf=gifti(LPGfp);
PG_LH=LPGf.cdata(:,1);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
RPGf=gifti(RPGfp);
PG_RH=RPGf.cdata(:,1);

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





% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=3;






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
% scale to 1
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
roybigbl_cm=roybigbl_cm(15:240,:);



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
for i=194:205
OpFlVecofInt=i;
TRofInt=validTRs(OpFlVecofInt);
u=OpFl.vf_right{OpFlVecofInt};
vATTR=fr.TRs{TRofInt};
% z-score
%vATTR=zscore(vATTR);
figure('units','pixels','position',[0 0 3500 3500])
axis([-1, 1, -1, 1, 0, 1]);
%quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), u(:, 1), u(:, 2), u(:, 3), 2, 'w');
%%%% souped up vectors
ret = bsxfun(@rdivide, u, sqrt(sum(u'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3)],'w',.7,'arrowRadius',.05)
%%%%%%
hold on
% for OpFl Vecs on PG
%trisurf(faces_r, vx_r(:, 1), vx_r(:, 2), vx_r(:, 3), PG_RH, 'EdgeColor','none');
%caxis([-5.5,6.5]);
% for OpFl Vecs on BOLD - divide by scaling factor?
%%% overlay pgg for ref angle clarity
PGG_ret=bsxfun(@rdivide, PGg_R, sqrt(sum(PGg_R'.^2))');
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[PGG_ret(g_noMW_combined_R,1), PGG_ret(g_noMW_combined_R,2), PGG_ret(g_noMW_combined_R,3)],'b',.7,'arrowRadius',.05)
hold on
quiver3D([Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor],[-PGG_ret(g_noMW_combined_R,1), -PGG_ret(g_noMW_combined_R,2), -PGG_ret(g_noMW_combined_R,3)],'r',.7,'arrowRadius',.05)
trisurf(faces_r, vx_r(:, 1)/scalingfactor, vx_r(:, 2)/scalingfactor, vx_r(:, 3)/scalingfactor, vATTR, 'EdgeColor','none');
caxis([-200,200])
axis equal
daspect([1, 1, 1]);
%colormap(roybigbl_cm);
colormap(roybigbl_cm);
c=colorbar;
c.FontSize=55;
c.LineWidth=3;
%c.Ticks=[-3 -2 -1 0 1 2 3];
c.Location='southoutside';
c.FontName='Arial';
view(280,185);;
%view(60,190)
i
fn=['~/boldvec_v2_' num2str(i) '_PGGoverlay.png'];
print(fn,'-dpng')
end

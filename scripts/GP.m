function GP(subj,imgFP,PGFP,TRstart,TRend,task,CorVal)
% read image
img=read_cifti(fp);

img=read_cifti('/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/sub-NDARINVUWH8GXXW/ses-baselineYear1Arm1/func/sub-NDARINVUWH8GXXW_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii')


% extract left hemi
CL_ind=img.diminfo{1}.models{1}.start+img.diminfo{1}.models{1}.vertlist;
CL=img.cdata(CL_ind,:);
% and right
CR_ind=img.diminfo{1}.models{2}.start+img.diminfo{1}.models{2}.vertlist;
CR=img.cdata(CR_ind,:);
CFull=[CL;CR];
% load in PG and arrange according to PG
pg=read_cifti(PGFP);
% extract left hemi
PGL_ind=pg.diminfo{1}.models{1}.start+pg.diminfo{1}.models{1}.vertlist;
PGL=pg.cdata(PGL_ind,:);
% and right
PGR_ind=pg.diminfo{1}.models{2}.start+pg.diminfo{1}.models{2}.vertlist;
PGR=pg.cdata(PGR_ind,:);
PGFull=[PGL;PGR];
% organize vertices in terms of their position on the PG
[Sorted,PGindices]=sort(PGFull(:,1));
% sort CFull
CFull=CFull(PGindices,:);


% L and R sep.
CL_ind=img.diminfo{1}.models{1}.start+img.diminfo{1}.models{1}.vertlist;
CL=img.cdata(CL_ind,:);
% and right
CR_ind=img.diminfo{1}.models{2}.start+img.diminfo{1}.models{2}.vertlist;
CR=img.cdata(CR_ind,:);
% load in PG and arrange according to PG
pg=read_cifti(PGFP);
% extract left hemi
PGL_ind=pg.diminfo{1}.models{1}.start+pg.diminfo{1}.models{1}.vertlist;
PGL=pg.cdata(PGL_ind,:);
% and right
PGR_ind=pg.diminfo{1}.models{2}.start+pg.diminfo{1}.models{2}.vertlist;
PGR=pg.cdata(PGR_ind,:);
% organize vertices in terms of their position on the PG
[Sorted,PGindicesL]=sort(PGL(:,1));
[Sorted,PGindicesR]=sort(PGR(:,1));
% sort CFull
CL=CL(PGindicesL,:);
CR=CR(PGindicesR,:);
% each vertex normalized w/r/t global mean and SD
% L
Lstds=mean(std(CL,0,2));
CL=CL-mean(CL,2);
CL=CL./Lstds;
% R
Rstds=mean(std(CR,0,2));
CR=CR-mean(CR,2);
CR=CR./Rstds;
filename=['/cbica/projects/abcdfnets/results/wave_output/' subj '_' task '_GP_' CorVal '.png'];
limz=[-3 3];


f = figure('visible','off');
subplot(2,1,1);
f=imagesc(CL,limz);


subplot(2,1,2);
f=imagesc(CR,limz);

colormap(gray);
%%%%%
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 15]); %x_width=10cm y_width=15cm
saveas(f,filename)







Cstds=mean(std(CFull,0,2));
CCFull=CFull-mean(CFull,2);
CFull=CFull./stds;


filename=['/cbica/projects/abcdfnets/results/wave_output/' subj '_' task '_GP_' CorVal '.png'];
f = figure('visible','off');
limz=[-4 4];
f=imagesc(CFull,limz);
colormap(gray);
colorbar








% each vertex normalized individually
%stds=std(CFull,0,2);
%CFull=CFull-mean(CFull,2);
%CFull=CFull./stds;

% each vertex normalized w/r/t global mean and SD
stds=mean(std(CFull,0,2));
CFull=CFull-mean(CFull,2);
CFull=CFull./stds;

filename=['/cbica/projects/abcdfnets/results/wave_output/' subj '_' task '_GP_' CorVal '.png'];
f = figure('visible','off');
limz=[-4 4];
f=imagesc(CFull,limz);
colormap(gray);
colorbar


%%%%%
set(f, 'PaperUnits', 'centimeters');
set(f, 'PaperPosition', [0 0 10 35]); %x_width=10cm y_width=15cm
saveas(f,filename)

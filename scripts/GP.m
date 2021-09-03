function GP(subj,imgFP,PGFP,TRstart,TRend,task,CorVal)
% read image
img=read_cifti(fp);
% extract left hemi
CL_ind=img.diminfo{1}.models{1}.start+img.diminfo{1}.models{1}.vertlist;
CL=img.cdata(CL_ind,TRstart:TRend);
% and right
CR_ind=img.diminfo{1}.models{2}.start+img.diminfo{1}.models{2}.vertlist;
CR=img.cdata(CR_ind,TRstat:TRend);
CFull=[CL;CR];
% load in PG and arrange according to PG
pg=read_cifti(PGFP);
% extract left hemi
PGL_ind=pg.diminfo{1}.models{1}.start+pg.diminfo{1}.models{1}.vertlist;
PGL=pg.cdata(PGL_ind,TRstart:TRend);
% and right
PGR_ind=pg.diminfo{1}.models{2}.start+pg.diminfo{1}.models{2}.vertlist;
PGR=pg.cdata(PGR_ind,TRstat:TRend);
PGFull=[PGL;PGR];
% organize vertices in terms of their position on the PG
[Sorted,PGindices]=sort(PGFull(:,1));
% sort CFull
CFull=CFull(PGindices,:);
stds=std(CFull,0,2);
CFull=CFull-mean(CFull,2);
CFull=CFull./stds;
filename=['/cbica/projects/abcdfnets/results/wave_output/' subj '_' task '_GP_' CorVal '.png'];
f = figure('visible','off');
limz=[-3 3];
f=imagesc(CFull,limz);
colormap(gray);
colorbar
saveas(f,filename)

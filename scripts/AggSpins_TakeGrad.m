% load i spins, iteratively project to surface and take gradient, save all to .mat

% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load in spins
sp_fp='/cbica/projects/pinesParcels/results/aggregated_data/PGPermuts_fs4.mat';
sp=load(sp_fp)

% load in surface data
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
% extract xyz coords for fs4, convert to spherical coords
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% initialize output file: elevation and azimuth directions, for each spin, for each face
spGPGGs=zeros(2,1000,10240);

% for s in 1000
for s=1:1000
	s
	% Extract spun elements
	sp_gPG_LH=sp.bigrotl(s,:);
	sp_gPG_RH=sp.bigrotr(s,:);
	% medial wall to 0
	sp_gPG_LH(sp_gPG_LH==100)=0;
	sp_gPG_RH(sp_gPG_RH==100)=0;
	% calculate spun group PG gradient on sphere
	sp_gPGg_L = grad(F_L, V_L, sp_gPG_LH);
	sp_gPGg_R = grad(F_R, V_R, sp_gPG_RH);
	% convert to elevation and azimuth to reduce output file size for iterative loading
	gazes_L=zeros(1,length(azd_L));
	gels_L=zeros(1,length(eld_L));
	for i=1:length(azd_L)
	    gvs_L=cart2sphvec(double([sp_gPGg_L(i,1);sp_gPGg_L(i,2);sp_gPGg_L(i,3)]),azd_L(i),eld_L(i));
	    gazes_L(i)=gvs_L(1);
	    gels_L(i)=gvs_L(2);
	    % drop the third vector, as each point is equidistant from the center of the sphere
	end
	% right hemi
	gazes_R=zeros(1,length(azd_R));
	gels_R=zeros(1,length(eld_R));
	for i=1:length(azd_R)
	    gvs_R=cart2sphvec(double([sp_gPGg_R(i,1);sp_gPGg_R(i,2);sp_gPGg_R(i,3)]),azd_R(i),eld_R(i));
	    gazes_R(i)=gvs_R(1);
	    gels_R(i)=gvs_R(2);
	end
	% insert left
	spGPGGs(1,s,1:5120)=gazes_L;
	spGPGGs(2,s,1:5120)=gels_L;
	% insert left
        spGPGGs(1,s,5121:10240)=gazes_R;
        spGPGGs(2,s,5121:10240)=gels_R;
end
% save output
save('/cbica/projects/pinesParcels/results/aggregated_data/PGGPermuts_fs4.mat','spGPGGs')

% plot to check
figure('units','pixels','position',[0 0 1000 1000])
axis([-1, 1, -1, 1, 0, 1]);
quiver3(P_L(:, 1), P_L(:, 2), P_L(:, 3), sp_gPGg_L(:,1), sp_gPGg_L(:,2), sp_gPGg_L(:,3), 4, 'k');
hold on
trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), sp_gPG_LH, 'EdgeColor','none');
axis equal
daspect([1, 1, 1]);
colorbar
view(180,60);
print('spunGPGG4_2.png','-dpng')

% needed paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
Subjs=readtable('~/PWs/rs_subs.csv')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load in surface data %%%%%
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4';
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
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

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

% load pgg
gpg=load('~/data/gpg_fs4.mat');
gPG_LH=gpg.gpg.gPG_LH;
gPG_RH=gpg.gpg.gPG_RH;
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% cartesian components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% true gazes and gels
tgazes_L=zeros(1,length(azd_L));
tgels_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    gvs_L=cart2sphvec(double([gPGx_L(i);gPGy_L(i);gPGz_L(i)]),azd_L(i),eld_L(i));
    tgazes_L(i)=gvs_L(1);
    tgels_L(i)=gvs_L(2);
	% drop the third vector, as each point is equidistant from the center of the sphere
end
% right hemi
tgazes_R=zeros(1,length(azd_R));
tgels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
    tgazes_R(i)=gvs_R(1);
    tgels_R(i)=gvs_R(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output array
subj_dip_ps=zeros(height(Subjs),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:height(Subjs)
	s
	Subj=table2array(Subjs(s,2));
	subj=Subj{:};
	% load subj data
	OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
	data=load(OpFlFp)
	% pull hemisphere vector fields
	vfl=data.us.vf_left;	
	vfr=data.us.vf_right;
	% get length of OpFl pairs
	lenOpFl=length(data.us.vf_left);
	% convert vectors
	azesOpf_L=zeros(length(azd_L),lenOpFl);
	elsOpf_L=zeros(length(azd_L),lenOpFl);
	for i=1:length(azd_L)
	    for fr=1:lenOpFl
		% current vector field
	        relVf_L=vfl{fr};
		% xyz components
	        xComp_L=relVf_L(i,1);
	        yComp_L=relVf_L(i,2);
	        zComp_L=relVf_L(i,3);
		% convert to spherical coord system
	        vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(i),eld_L(i));
	        % store in output vector (r is redundant across all vecs)
		azesOpf_L(i,fr)=vs_L(1);
	        elsOpf_L(i,fr)=vs_L(2);
	        rvec=vs_L(3);
	    end
	end
	% right hemisphre
	azesOpf_R=zeros(length(azd_R),lenOpFl);
	elsOpf_R=zeros(length(azd_R),lenOpFl);
	for i=1:length(azd_R)
	    for fr=1:lenOpFl
		% current vector field
	        relVf_R=vfr{fr};
		% xyz components
	        xComp_R=relVf_R(i,1);
	        yComp_R=relVf_R(i,2);
	        zComp_R=relVf_R(i,3);
		% convert to spherical coord system
	        vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(i),eld_R(i));
	        % store in output vector (r is redundant across all vecs)
		azesOpf_R(i,fr)=vs_R(1);
	        elsOpf_R(i,fr)=vs_R(2);
	        rvec=vs_R(3);
	    end
	end
	disp('done converting opfl vectors from cartesian')

	% calculate angular distances

	True_gangDist_L=zeros(lenOpFl,length(azd_L));
	True_gangDist_R=zeros(lenOpFl,length(azd_R));
	% for each vertex
	for Vert=1:length(azd_L)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_L=[tgazes_L(Vert) tgels_L(Vert)]; 
	    % PG GROUP LOAD IN
	    for fr=1:lenOpFl
	        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
	    	ga = acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )));
		True_gangDist_L(fr,Vert) = ga;	
	    end
	end
	for Vert=1:length(azd_R)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_R=[tgazes_R(Vert) tgels_R(Vert)];
	    % PG GROUP LOAD IN
	    for fr=1:lenOpFl
	        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
	        ga = acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )));
	        True_gangDist_R(fr,Vert) = ga;
	    end
	end

	% make the last one "true" observed from masked data
	Tr_AngDists=horzcat(True_gangDist_L(:,g_noMW_combined_L),True_gangDist_R(:,g_noMW_combined_R));
	[dip, p_value, xlow,xup]=HartigansDipSigniftest(Tr_AngDists,10000)
	subj_dip_ps(s)=p_value;
end
writetable(table(subj_dip_ps),'~/results/PWs/native_dip_ps.csv')


function PGG_AngDistCalc4_sp(subj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load in optical flow data and calculate gradient gradient data, compare angular distance between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in fsav4 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
data=load(OpFlFp)

%%% Load in surface data
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
% vector fields
vfl=data.us.vf_left;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% vector fields
vfr=data.us.vf_right;
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

% load in spun PGGs
sp_PGGs=load('/cbica/projects/pinesParcels/results/aggregated_data/PGGPermuts_fs4.mat');
gazes_L_all=squeeze(sp_PGGs.spGPGGs(1,:,1:5120));
gels_L_all=squeeze(sp_PGGs.spGPGGs(2,:,1:5120));
gazes_R_all=squeeze(sp_PGGs.spGPGGs(1,:,5121:10240));
gels_R_all=squeeze(sp_PGGs.spGPGGs(2,:,5121:10240));

% get length of OpFl pairs
lenOpFl=length(data.us.vf_left);

% translate xyz vector fields from opfl to az/el/r
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

% actually calculate Opflow vector vs PGG vector angular distance
gangDist_L=zeros(lenOpFl,length(azd_L),1000);
gangDist_R=zeros(lenOpFl,length(azd_R),1000);

% for each spin
for S=1:1000
	tic
	S
	% this iterations gazels and gel
	gazes_L=gazes_L_all(S,:);
	gels_L=gels_L_all(S,:);
	gazes_R=gazes_R_all(S,:);
	gels_R=gels_R_all(S,:);
	% for each vertex
	for Vert=1:length(azd_L)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_L=[gazes_L(Vert) gels_L(Vert)]; 
	    % PG GROUP LOAD IN
	    for fr=1:lenOpFl
	        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
    		ga = acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )));
		gangDist_L(fr,Vert,S) = ga;	
    	     end
	end
	toc
	tic
	% right hemi
	% for each vertex
	for Vert=1:length(azd_R)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_R=[gazes_R(Vert) gels_R(Vert)];
	    for fr=1:lenOpFl
	        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
        	ga = acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )));
        	gangDist_R(fr,Vert,S) = ga;
    	    end
	end
	toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this needs to be dip tested in matlab. 1500*5120*1000 output files will not jive well with storage constraints
% which means it needs to be MASKED first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save it out
%AngDist=struct;
%AngDist.gLeft=gangDist_L;
%AngDist.gRight=gangDist_R;
%AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat4_sp.mat'];
%save(AngDistFP,'AngDist')

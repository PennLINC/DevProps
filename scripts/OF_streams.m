function OF_streams(subj)
%addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% enter subj manually for now
%subj='';

%%% load in OF vector fields
OFvecsFile=load(['~/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat']);
OFvecs=OFvecsFile.us;
OFvecs_L=OFvecs.vf_left;
OFvecs_R=OFvecs.vf_right;


%%% load in spherical coordinates
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
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

%%% create cartesian sphere grid for streamline
% density of grid = n
n=40;
% interpolate vector field to meshgrid:
% neighboring grid points should be 9 units away from each other, 
MaxCoordOfSph=110;
% we want points only on the surface, so threshold anything that's not pretty dang close to outer "shell"
MinCoordOfSph=90;
% create initial cartesian meshgrid:
%https://de.mathworks.com/matlabcentral/answers/498723-how-to-start-streamlines-on-the-surface-of-a-sphere-griddedinterpolant-requires-at-least-two-sample
x = linspace((-MaxCoordOfSph+10),(MaxCoordOfSph-10),n);
[X,Y,Z] = meshgrid(x,x,x);
% field
Fx = X;
Fy = Y;
Fz = Z;
% remove inside of sphere
idx = sqrt(X.^2+Y.^2+Z.^2) < MinCoordOfSph;
Fx(idx) = NaN;
Fy(idx) = NaN;
Fz(idx) = NaN;
% remove outside of sphere
idx = sqrt(X.^2+Y.^2+Z.^2) > MaxCoordOfSph;
Fx(idx) = NaN;
Fy(idx) = NaN;
Fz(idx) = NaN;

%%% pull in mask
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

%%% initialize loop variables
% adjacency matrix for streamlines
adjMatL=zeros(5120);
adjMatR=zeros(5120);

% sequence of numbers for indexing faces
faceSeq=1:5120';

% for all TR pairs
for trp=1:length(OFvecs_L)
tic
trp
%%%%%%%%%%%%%
%%%%%%%% LEFT
%%%%%%%%%%%%%

	%%% select vector field to propagate along
	VecField=OFvecs_L{trp};

	% using griddata for 3d
	uq = griddata(P(:,1),P(:,2),P(:,3),VecField(:,1),Fx,Fy,Fz);
	vq = griddata(P(:,1),P(:,2),P(:,3),VecField(:,2),Fx,Fy,Fz);
	wq = griddata(P(:,1),P(:,2),P(:,3),VecField(:,3),Fx,Fy,Fz);

	% replace NaNs with 0s
	uq(isnan(uq))=0;
	vq(isnan(vq))=0;
	wq(isnan(wq))=0;

	% generate streamlines using only "in-bounds" faces as seed locations
	trp_streamlines=stream3(X,Y,Z,uq,vq,wq,P(g_noMW_combined_L,1),P(g_noMW_combined_L,2),P(g_noMW_combined_L,3),[.1,500]);

	% for each streamline, find closest point. Reduce to sequence/vector of points falling along this streamline (instead of xyz coords)
	for f=1:length(g_noMW_combined_L)
		% what is true index of this f
		fInd=g_noMW_combined_L(f);
		% recover potential streamlines from this face
		Streams=trp_streamlines{:,f};
		sizeStreams=size(Streams);
		% initialize "closest" vector to record faces the streamlines grace
		closestVec=zeros(1,sizeStreams(1));
		% for each segment of stream (reminder that 1st is seed point)
		for segm=1:sizeStreams(1);
			% get closest neighbors to the described points on freesurfer mesh centroids
			% https://www.mathworks.com/matlabcentral/answers/107029-find-closest-coordinates-to-a-point
			dist2 = sum((P - Streams(segm,:)) .^ 2, 2);
			closest = faceSeq(dist2 == min(dist2));
			% -1 to account for first point being "lost" (streamline in the seed-point itself)
			closestVec(segm)=closest;
		end
		% get unique faces in closest vector
		unFaces=unique(closestVec);
		% get non-self values
		nonSelf=setdiff(unFaces,fInd);
		% if there are non-self values, insert into adjacency
		if length(nonSelf) > 0
			adjMatL(fInd,nonSelf)=adjMatL(fInd,nonSelf)+1;
		end
	end

%%%%%%%%%%%%%
%%%%%%% RIGHT
%%%%%%%%%%%%%


        %%% select vector field to propagate along
        VecField=OFvecs_R{trp};

	% using griddata for 3d
        uq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),VecField(:,1),Fx,Fy,Fz);
        vq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),VecField(:,2),Fx,Fy,Fz);
        wq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),VecField(:,3),Fx,Fy,Fz);

        % replace NaNs with 0s
        uq(isnan(uq))=0;
        vq(isnan(vq))=0;
        wq(isnan(wq))=0;

        % generate streamlines using only "in-bounds" faces as seed locations
        trp_streamlines=stream3(X,Y,Z,uq,vq,wq,Pr(g_noMW_combined_R,1),Pr(g_noMW_combined_R,2),Pr(g_noMW_combined_R,3),[.1,500]);
	
	% for each streamline, find closest point. Reduce to sequence/vector of points falling along this streamline (instead of xyz coords)
        for f=1:length(g_noMW_combined_R)
		% what is true index of this f
                fInd=g_noMW_combined_R(f);
                % recover potential streamlines from this face
                Streams=trp_streamlines{:,f};
                sizeStreams=size(Streams);
                % initialize "closest" vector to record faces the streamlines grace
                closestVec=zeros(1,sizeStreams(1));
                % for each segment of stream (reminder that 1st is seed point)
                for segm=1:sizeStreams(1);
                        % get closest neighbors to the described points on freesurfer mesh centroids
                        % https://www.mathworks.com/matlabcentral/answers/107029-find-closest-coordinates-to-a-point
                        dist2 = sum((Pr - Streams(segm,:)) .^ 2, 2);
                        closest = faceSeq(dist2 == min(dist2));
                        % -1 to account for first point being "lost" (streamline in the seed-point itself)
                        closestVec(segm)=closest;
                end
                % get unique faces in closest vector
                unFaces=unique(closestVec);
                % get non-self values
                nonSelf=setdiff(unFaces,fInd);
                % if there are non-self values, insert into adjacency
                if length(nonSelf) > 0
                        adjMatR(fInd,nonSelf)=adjMatR(fInd,nonSelf)+1;
                end
        end

% end over this TRP %%%%
toc %%%%%%%%%%%%%%%%%%%
end %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% mask mw
adjMatL=adjMatL(g_noMW_combined_L,g_noMW_combined_L);
adjMatR=adjMatR(g_noMW_combined_R,g_noMW_combined_R);

% saveout
adjMats=struct('L',{adjMatL},'R',{adjMatR});
fn=['~/results/PWs/Proced/' subj '/' subj '_FSC.mat'];
save(fn,'adjMats')
 

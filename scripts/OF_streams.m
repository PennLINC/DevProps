addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% enter subj manually for now
subj='';

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

%%%% now we need to convert vector fields to this space

%%% select example vector field to propagate along
exVecField=OFvecs_R{1059};

% using griddata for 3d
uq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),exVecField(:,1),Fx,Fy,Fz);
vq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),exVecField(:,2),Fx,Fy,Fz);
wq = griddata(Pr(:,1),Pr(:,2),Pr(:,3),exVecField(:,3),Fx,Fy,Fz);

% replace NaNs with 0s
uq(isnan(uq))=0;
vq(isnan(vq))=0;
wq(isnan(wq))=0;

% viz this versus original
%figure
%quiver3(Fx,Fy,Fz,uq,vq,wq,2,'r');
%print('~/regrid_quiv.png','-dpng');
%figure
%quiver3(Pr(:, 1), Pr(:, 2), Pr(:, 3), exVecField(:, 1), exVecField(:, 2), exVecField(:, 3), 2, 'b');
%print('~/OG_quiv.png','-dpng');

%%% calculate streamlines
% x y z are COORDINATES of vector field (incenters of triangle locations in cartesian grid)
% U V W are the vector fields themselves
% to seed from every face, startX,startY,startZ should be same as x y z 

figure('units','pixels','position',[0 0 1500 1500])
examplestreamlines=stream3(X,Y,Z,uq,vq,wq,Pr(:,1),Pr(:,2),Pr(:,3),[.1,500]);
streamline(examplestreamlines)

% sequence of numbers for indexing faces
faceSeq=1:5120;

% for each streamline, find closest point. Reduce to sequence/vector of points falling along this streamline (instead of xyz coords)
for f=1:5120
	% recover potential streamlines from this face
	Streams=examplestreamlines{:,f};
	% for each segment of stream (1st is seed point)
	for segm=2:length(Streams)
		% get closest neighbors to the described points on freesurfer mesh centroids
		dist2 = sum((Pr - Streams(segm,:)) .^ 2, 2);
		closest = Pr(dist2 == min(dist2),:);
	end
f
end

view(280,185)
print('~/test3.png','-dpng')


% for each streamline point, find closest point
% https://www.mathworks.com/matlabcentral/answers/107029-find-closest-coordinates-to-a-point
% record instances in "SC' matrix

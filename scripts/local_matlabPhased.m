data=load('~/data.mat')
addpath(genpath('~/ofd'))

% faces_L
F=data.data.faces_l;
% vertices V
V=data.data.vx_l;
% vector fields
vfl=data.data.vf_left;

% get incenters of triangles
TR = TriRep(F,V);
P = TR.incenters;

% load in PG
PG_LH=load('~/PG_LH.mat')
PG_LH=PG_LH.PG_LH;
% plot PG_LH
trisurf(F, V(:, 1), V(:, 2), V(:, 3), PG_LH)

% calculate PG gradient on sphere
PGg = grad(F, V, PG_LH);

% extract face-wise vector cartesian vector components
PGx=PGg(:,1);
PGy=PGg(:,2);
PGz=PGg(:,3);

% plot gradient gradient vector field
trisurf(F, V(:, 1), V(:, 2), V(:, 3), PG_LH,'EdgeColor', 'none')
hold on
quiver3(P(:, 1), P(:, 2), P(:, 3), PGx, PGy, PGz,'b');

% convert angles to format for cart2sphvec
Cart_angleMat=vertcat(PGx',PGy',PGz');

% translate xyz spherical coordinates to az/el/r
[az,el,r]=cart2sph(P(:,1),P(:,2),P(:,3));

% convert to degrees from radians
azd=rad2deg(az);
eld=rad2deg(el);

% translate xyz vector components at coordinates to az/el/r
azes=[];
els=[];
for i=1:length(azd)
    vs=cart2sphvec(double([PGx(i),PGy(i),PGz(i)]'),azd(i),eld(i));
    azes(i)=vs(1);
    els(i)=vs(2);
    rvec=vs(3)
end

% translate xyz vector fields from opfl to az/el/r
azesOpf=zeros(length(azd),1811);
elsOpf=zeros(length(azd),1811);
thetasOpf=zeros(length(azd),1811);
for i=1:length(azd)
    i
    for fr=1:1811
        relVf=vfl{fr};
        xComp=relVf(i,1);
        yComp=relVf(i,2);
        zComp=relVf(i,3);
        vs=cart2sphvec(double([xComp,yComp,zComp]'),azd(i),eld(i));
        azesOpf(i,fr)=vs(1);
        elsOpf(i,fr)=vs(2);
        rvec=vs(3);
    end
end

% convert pseudo y x to single angle (theta), az as x, el as y
PGg_thetas=atan2(els,azes);
% loop over each vertex
for i=1:length(azd)
    Opfl_thetas(i,:)=atan2(elsOpf(i,:),azesOpf(i,:));
end

% get mode and compare to PG

% initialize BU array
BU_angDist=zeros(1811,length(azd));

% for each vertex
for Vert=1:length(azd)
    PGvec=[azes(Vert) els(Vert)];
    OpFlVecs=horzcat(azesOpf(Vert,:)',elsOpf(Vert,:)');
    for fr=1:1811
        BUCosTheta = max(min(dot(PGvec,OpFlVecs(fr,:))/(norm(PGvec)*norm(OpFlVecs(fr,:))),1),-1);
        BU_angDist(fr,Vert) = real(acosd(BUCosTheta));
    end
end


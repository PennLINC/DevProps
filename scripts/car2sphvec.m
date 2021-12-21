%   translates a vector in rectangular coordinates to its corresponding
%   spherical representation in the coordinates (az_hat, el_hat, r_hat )
%   defined in the direction of (az;el).

% Ar is 3xN rectangular coordinates [x;y;z]
% As is 3xN spherical coordinates [az;el;r] (see cart2sph)
% az is x->y in degrees, scalar
% el is xy->z in degrees, scalar

eml_assert_no_varsize(2:3,Ar,az,el);
validateattributes(Ar,{'double'},{'finite','nonnan','nonempty','2d','nrows',3},'cart2sphvec','VR');
sigdatatypes.validateAngle(az,'cart2sphvec','AZ',{'scalar','<=',180,'>=',-180});
sigdatatypes.validateAngle(el,'cart2sphvec','EL',{'scalar','<=',90,'>=',-90});
M = phased.internal.azelcoord(az,el);

As = M.'*Ar;


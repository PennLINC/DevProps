function V = VecNormalize(V)
%NORMALIZE Normalizes row vectors.
%
%   V = normalise(V) takes a n-by-3 matrix V and returns normalises row
%   vectors.
len = sqrt(sum(V.^2, 2));
V = bsxfun(@rdivide, V, len);
end

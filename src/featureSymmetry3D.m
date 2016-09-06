function [FS,FA] = featureSymmetry3D(m1,m2,m3,m4,T)
%
%     [FS,FA] = featureSymmetry3D(m1,m2,m3, m4, T = 0.18)
%
%  Calculates the feature symmetry (FS) and feature asymmetry (FA) using
%  the components of the monogenic signal (m1,m2,m3,m4). These arrays may be
%  four dimensional, with the fourth dimension corresponding to different
%  wavelengths of the bandpass filter. T is the threshold used.
%
%  Adapted by Chris Bridge (March 2014) from code by Vicente Grau and Ana
%  Namburete
%  christopher.bridge@eng.ox.ac.uk

% Threshold
if nargin < 5
    T = 0.18;
end

[ysize, xsize, zsize, ssize] = size(m1);

% Small constant to avoid division by zero
epsilon = 0.001;

% Combine the odd components
odd = sqrt(m2.^2 + m3.^2 + m4.^2);
even = abs(m1);

% Calculate the denominator (= local energy + epsilon)
denominator = sqrt(even.^2 + odd .^2) + epsilon;

% Calculate the numerators for FA and FS at all scales
% NB no need to take absolute of 'odd' as it must be positive due to the
% way it's calculated
FS_numerator = max(even - odd - T, zeros(ysize, xsize, zsize, ssize));
FA_numerator = max(odd - even - T, zeros(ysize, xsize, zsize, ssize));

% Divide numerator by denominator
FS = FS_numerator./denominator; 
FA = FA_numerator./denominator;

% Sum across scales, and divide by nimber of scales to give value between 0
% and 1 (i.e. take mean across the scale dimension)
FS = mean(FS, 4);
FA = mean(FA, 4);


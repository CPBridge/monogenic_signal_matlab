function [FA_y, FA_x, FA_z, FS] = orientedSymmetry3D(m1,m2,m3,m4,T)
%
%  [FA_y, FA_x, FA_z, FS] = orientedAsymmetry3D(m1,m2,m3,m4,T = 0.18)
%
%   Creates an edge map using oriented feature assymetry calculated from
%   the monogenic signal (m1,m2,m3,m4). FA_y contains the vertical edge
%   components and FA_x contains the horizontal edge components. The
%   direction of the edges is defined to point in the direction of
%   increasing intensity, using the standard MATLAB definitions for the
%   increasing y and x directions (i.e. down and right).
%   
%   FS returns the oriented feature symmetry map (i.e. feature symmetry
%   with a sign to indiciate the polarity of the symmetry)

% Threshold
if nargin < 5
    T = 0.18;
end

[ysize, xsize, tsize, ssize] = size(m1);

% Small constant to avoid division by zero
epsilon = 0.001;

% Take the absolute value of m1 as the even part of the filter
even = abs(m1);
odd = sqrt(m2.*m2 + m3.*m3 + m4.*m4);

% Calculate the denominator (= local energy + epsilon)
denominator = sqrt(even.*even + m2.*m2 + m3.*m3 + m4.*m4) + epsilon;

% Calculate the numerators for FA and FS at all scales
% NB no need to take absolute value of 'odd' as it must be positive due to
% the way it's calculated
FA_numerator = max(odd - even - T, zeros(ysize, xsize, tsize, ssize));
FS_numerator = max(even - odd - T, zeros(ysize, xsize, tsize, ssize));

% Divide numerator by denominator and include orientation
FA_y = (FA_numerator./denominator).*(m2./odd);
FA_x = (FA_numerator./denominator).*(m3./odd);
FA_z = (FA_numerator./denominator).*(m4./odd);
FS = (FS_numerator./denominator).*sign(m1);

% Sum across scales, and divide by nimber of scales to give value between 0
% and 1 (i.e. take mean across the scale dimension)
FA_y = mean(FA_y, 4);
FA_x = mean(FA_x, 4);
FA_z = mean(FA_z, 4);
FS = mean(FS, 4);


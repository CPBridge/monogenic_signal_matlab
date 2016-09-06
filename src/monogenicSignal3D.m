function [Fm1, Fm2, Fm3, Fm4] = monogenicSignal3D(vol, filtStruct)
%
%        [Fm1, Fm2, Fm3, Fm4] = monogenicSignal3D(vol, filtStruct)
%
% This function computes the monogenic signal using the Riesz transform as
% suggested by Felsberg.
%
%   Inputs:
%       vol        = volume for which monogenic signal will be computed
%       filtStruct = a structure containing the necessary filters, as
%                    returned by createMonogenicFilters3D
%
%   Outputs:
%       Fm1             = even component of the monogenic signal
%       Fm2,Fm3, Fm4    = odd components of the monogenic signal
%
% Adapted by Chris Bridge (March 2014) from code by Vincente Grau and Ana
% Namburete
% christopher.bridge@eng.ox.ac.uk

% Create output arrays
outputsize = cat(2, size(vol), filtStruct.numFilt);
Fm1 = zeros(outputsize);
Fm2 = zeros(outputsize);
Fm3 = zeros(outputsize);
Fm4 = zeros(outputsize);

% Compute the 3-dimensional fast Fourier transform of the original image
F = fftn(vol);

% Filter using the Reisz filter
R_03 = F.*filtStruct.ReiszFilt03;
R_12 = F.*filtStruct.ReiszFilt12;

% Compute the parts of the monogenic signal
for flt = 1:filtStruct.numFilt
    F03 = ifftn(R_03.*filtStruct.bpFilt(:,:,:,flt));
    F12 = ifftn(R_12.*filtStruct.bpFilt(:,:,:,flt));
    Fm1(:,:,:,flt) = real( F03 );
    Fm2(:,:,:,flt) = real( F12 );
    Fm3(:,:,:,flt) = imag( F12 );
    Fm4(:,:,:,flt) = imag( F03 );
end

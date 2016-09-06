function [Fm1, Fm2, Fm3] = monogenicSignal(im, filtStruct)
%
%        [Fm1, Fm2, Fm3] = monogenicSignal(im, filtStruct)
%
% This function computes the monogenic signal using the Riesz transform as
% suggested by Felsberg. This function returns the 2D monogenic signal of
% each frame regardless of the input dimensionality
%
%   Inputs:
%       im         = image for which monogenic signal will be computed.
%                    May be 2 or 3 dimensional stack of images.
%       filtStruct = a structure containing the necessary filters, as
%                    returned by createMonogenicFilters
%
%   Outputs:
%       Fm1        = even component of the monogenic signal
%       Fm2,Fm3    = odd components of the monogenic signal
%                    The first three dimensions of the output are the image
%                    dimensions, the fourth is the bandpass wavelength if
%                    multiple are used.
%
%
% Adapted by Chris Bridge (March 2014) from code by Vincente Grau and Ana
% Namburete
% christopher.bridge@eng.ox.ac.uk

% Compute the 2-dimensional fast Fourier transform of the original image or
% stack of images
F = fft2(im);
Ffilt = bsxfun(@times, F, filtStruct.bpFilt);

% Compute the parts of the monogenic signal (NB applying the ifft2 to a 3D
% image here performs the 2D ifft to each two dimensional slice)
Fm1 = real(ifft2( Ffilt ));    %even component
Fmodd = ifft2( bsxfun(@times, Ffilt, filtStruct.ReiszFilt) );
Fm2 = real(Fmodd);  %odd components...
Fm3 = imag(Fmodd);

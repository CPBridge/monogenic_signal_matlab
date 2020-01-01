function filtStruct = createMonogenicFilters3D(ysize,xsize,zsize,wl_xy,wl_z,filtType,parameter)
%
%   filtStruct = createMonogenicFilters3D(ysize, xsize, zsize , wl_xy, wl_z, filtType = 'lg', parameter)
%
%    Creates a set of frequency domain filters combining bandpass and Reisz
%    filters for use in calculating the monogenic signal.
%
%         Inputs:
%                ysize, xsize, zsize - Y and X dimensions of the filters
%                wl_xy    - vector of spatial wavelengths to use (for
%                           filtering image dimensions 1 and 2)
%                wl_z     - vector of temporal wavelengths to use (for
%                           filtering image dimension 3)
%                filtType - Filter type to use:
%                           'lg'    - Log-Gabor radial filter
%                           'gabor' - Gabor radial filter
%                           'gd'    - Gaussian derivative
%                           'cau'   - Cauchy
%                           'dop'   - difference of Poisson filter
%                 parameter - The shape parameters for 'lg' and 'dop' type
%                             filters
%
%         Outputs:
%                 fltStruct - Output structure containing the filters used
%                           for use in calculating the mongenic signal.
%
%                          Fields:
%                           -  bpFilt - A four-dimensional array. Stacked
%                              along the fourth dimension are the bandpass
%                              filters for each wavelength.
%                           -  ReiszFilt03, ReiszFilt12, two complex-valued Reisz
%                              transform filters that allow the monogenic signal
%                              to be computed with two IFFTs. The first is for the
%                              even and third (z) odd part, the second for the
%                              first and second odd parts (x,y)
%                           -  numFilt - The number of bandpass filters
%
%
% Adapted by Chris Bridge (March 2014) from code by Vicente Grau and Ana
% Namburete
% christopher.bridge@eng.ox.ac.uk

if nargin < 6
    filtType = 'lg';
end

% Parameters governing the bandwidth of the bandpass filter for Gabor and
% log-Gabor type filters, and ratio of centre frequencies for difference of
% Poisson filters
if nargin > 6
    if(parameter < 0.0 || parameter > 1.0)
        error('Parameter must be between 0.0 and 1.0');
    end
    sigmaOnf = parameter;
    ratio = parameter;
else
    sigmaOnf = 0.5;
    ratio = 0.98;
end

% Check that we have the same number of wavelength values for the spatial
% and temporal filters. If one is single valued, extend it to use the same
% value for each
if length(wl_xy) ~= length(wl_z)
    if length(wl_xy) == 1
        wl_xy = wl_xy*ones(length(wl_z));
    elseif length(wl_z) == 1
        wl_z = wl_z*ones(length(wl_xy));
    else
        error('The dimensions of the spatial and temporal wavelength vectors must agree, or one must be scalar')
    end
end

% Determine the number of filters to try from input
numFilt = length(wl_xy);

[yGrid, xGrid, zGrid] = freqgrid3(ysize,xsize,zsize);

% Output structure
bpFilt = zeros(ysize,xsize,zsize,numFilt);

% Generate band-pass filters (in frequency domain)
for flt = 1:numFilt
    % Construct the filter - first calculate the radial filter component.
    fo = 1.0/wl_xy(flt);   % Centre frequency of spatial filter.
    w0 = fo;       % Normalised radius from centre of frequency plane
                       % corresponding to fo.

    g0 = 1.0/wl_z(flt);   % Centre frequency of temporal filter.
    u0 = g0/0.5;       % Normalised radius from centre of frequency plane
                       % corresponding to fo.

    % Determine the spatial regions to use
    w = sqrt( (yGrid.^2)./(w0^2) + (xGrid.^2)./(w0^2) + (zGrid.^2)./(u0^2) );
    w(1,1,1) = 1;    %Trick to avoid division by zero

    if strcmp(filtType, 'lg')
        bpFilt(:,:,:,flt) = exp((-(log(w)).^2) / (2 * log(sigmaOnf)^2));         % computation of 3D log-gabor filter across the volumetric range
    elseif strcmp(filtType, 'gabor')
        bpFilt(:,:,:,flt) = exp((-(w).^2) / (2 * sigmaOnf^2));                   % gabor filter
    elseif strcmp(filtType, 'gd')
        bpFilt(:,:,:,flt) = w .* exp(-4*(w.^2));                          % isotropic gaussian derivative filter; wl is used for sigma
    elseif strcmp(filtType, 'cau')
        bpFilt(:,:,:,flt) = w .* exp(-2*w); % wl is used for sigma
    elseif strcmp(filtType, 'dop')
        s2 = log(ratio)/(ratio-1);
        s1 = ratio*s2;
        bpFilt(:,:,:,flt) = exp(-2*w*s1) - exp(-2*w*s2);                  % difference of Poisson filters
    end

    % Set the DC value of the filter
    bpFilt(1, 1, 1, flt) = 0;

    % Also remove unwanted high frequency components in filters
    % with even dimensions
    if(mod(ysize,2) == 0)
       bpFilt( ysize/2 +1,:, 1, flt) = 0;
    end
    if(mod(xsize,2) == 0)
       bpFilt( :, xsize/2 +1, 1, flt) = 0;
    end
    if(mod(zsize,2) == 0)
       bpFilt( :, :, zsize/2 +1, flt) = 0;
    end
end

% Normalise by the maximum value of the sum of all filters
sumFilt = sum(bpFilt, 4);
filtStruct.bpFilt = bpFilt ./ max(sumFilt(:));

% Generate the Riesz filter components (i.e. the odd filter whose
% components are imaginary)
w = sqrt(yGrid.^2 + xGrid.^2 + zGrid.^2);
w(1, 1, 1) = 1;
filtStruct.ReiszFilt03 = 1 - (zGrid ./ w);        % Complex filter for the even and z-dimension odd component
filtStruct.ReiszFilt12 = (1i*yGrid - xGrid)./ w; % Complex filter for the spatial odd components
filtStruct.numFilt = numFilt;

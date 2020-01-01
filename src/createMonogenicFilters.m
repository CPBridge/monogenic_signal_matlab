function filtStruct = createMonogenicFilters(ysize,xsize,wl,filtType,parameter)
%
%   filtStruct = createMonogenicFilters(ysize,xsize,wl, filtType = 'lg', parameter)
%
%    Creates a set of frequency domain filters combining bandpass and Reisz
%    filters for use in calculating the monogenic signal.
%
%         Inputs:
%                ysize, xsize - Y and X dimensions of the filters
%                wl           - Vector of wavelengths to use
%                filtType     - Filter type to use:
%                               'lg'    - Log-Gabor radial filter
%                               'gabor' - Gabor radial filter
%                               'gd'    - Gaussian derivative filter
%                               'cau'   - Cauchy
%                               'dop'   - Difference of Poisson filter
%                parameter    - The shape parameters for 'lg' and 'dop'
%                               type filters
%
%         Outputs:
%                 fltStruct   - Output structure containing the filters used
%                               for use in calculating the mongenic signal.
%
%                               Fields:
%                             -  bpFilt - A four-dimensional array. Stacked
%                              along the fourth dimension are the bandpass
%                              filters for each wavelength. The third
%                              dimension is singleton. This allows
%                              intuitive filtering of stacks of images
%                              (i.e. videos)
%                             -  ReiszFilt - the complex-valued filter for the
%                              the Reisz transform
%                             -  numFilt - The number of bandpass filters
%
%
% Adapted by Chris Bridge (March 2014) from code by Vicente Grau and Ana
% Namburete
% christopher.bridge@eng.ox.ac.uk

if nargin < 4 || isempty(filtType)
    filtType = 'lg';
end

% Frequency grid for the filter
[yGrid,xGrid] = freqgrid2(ysize,xsize);

% Determine the spatial regions to use
w = sqrt(yGrid.^2 + xGrid.^2);
w(1, 1) = 1;    %Trick to avoid division by zero at DC component

% Determine the number of filters to try from input
numFilt = length(wl);

% Parameters governing the bandwidth of the bandpass filter for Gabor and
% log-Gabor type filters, and ratio of centre frequencies for difference of
% Poisson filters

if nargin > 4 && ~isempty(parameter)
    if(parameter < 0.0 || parameter > 1.0)
        error('Parameter must be between 0.0 and 1.0');
    end
    sigmaOnf = parameter;
    ratio = parameter;
else
    sigmaOnf = 0.5;
    ratio = 0.98;
end

% Output structure
bpFilt = zeros(ysize,xsize,1,numFilt);


% Generate band-pass filters (in frequency domain)
for flt = 1:numFilt
    % Construct the filter - first calculate the radial filter component.
    fo = 1.0/wl(flt);   % Centre frequency of filter.
    w0 = fo;       % Normalised radius from centre of frequency plane
                        % corresponding to fo.
    if strcmp(filtType, 'lg')
        bpFilt(:,:,1,flt) = exp((-(log(w/w0)).^2) / (2 * log(sigmaOnf)^2));         % computation of 3D log-gabor filter across the volumetric range
    elseif strcmp(filtType, 'gabor')
        bpFilt(:,:,1,flt) = exp((-(w/w0).^2) / (2 * sigmaOnf^2));                   % gabor filter
    elseif strcmp(filtType, 'gd')
        bpFilt(:,:,1,flt) = w .* exp(-(w.^2)*(wl(flt)^2));                          % isotropic gaussian derivative filter; wl is used for sigma
    elseif strcmp(filtType, 'cau')
        bpFilt(:,:,1,flt) = w .* exp(-(w)*(wl(flt)));                               % cauchy filter, wl is used for sigma
    elseif strcmp(filtType, 'dop')
        s2 = wl(flt)/((ratio-1))*log(ratio);                                   % difference of Poisson filters
        s1 = ratio*s2;
        bpFilt(:,:,1,flt) = exp(-w*s1) - exp(-w*s2);
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
end

% Normalise by the maximum value of the sum of all filters
sumFilt = sum(bpFilt, 4);
filtStruct.bpFilt = bpFilt ./ max(sumFilt(:));
filtStruct.numFilt = numFilt;

% Generating the Riesz filter components as a complex filter
filtStruct.ReiszFilt = (1i * yGrid - xGrid)./ w;     % (iY) + i(iX) = iY - X

% Create a frequency-domain differentiation filter (using the
% differentation property of Fourier Transforms)
filtStruct.diffFilt = 2i*pi*yGrid - 2*pi*xGrid;          % (iY) + i(iX) = iY - X

filtStruct.wl = wl;
filtStruct.filtType = filtType;
filtStruct.sigmaOnf = sigmaOnf;
filtStruct.ratio = ratio;

function [yGrid, xGrid] = freqgrid2(ysize,xsize,gpu)

%  [yGrid, xGrid] = freqgrid2(ysize,xsize)
%  [yGrid, xGrid] = freqgrid2(size)
%
% Creates a 2D grid of frequencies of a given size, respecting the MATLAB
% convention for the location of the DC component (top left). This means
% that if a filter is desgined on this grid, ifft2(filter) will behave as
% expected. Use fftshift before visualising.
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

if nargin < 2
    xsize = ysize;
end

% Dimensions
ymid = floor(ysize/2); 
xmid = floor(xsize/2);

% Work out the maximum frequency in the grid - depends on whether the image
% has odd or even dimensions due to the definition of the FFT
if(mod(ysize,2)) == 0
    ymax = ymid-1;
else
    ymax = ymid;
end

if(mod(xsize,2)) == 0
    xmax = xmid-1;
else
    xmax = xmid;
end

% Create the grid for the filter
[yGrid, xGrid] = ndgrid( -ymid : ymax, -xmid : xmax);
yGrid = ifftshift(yGrid)/ysize;
xGrid = ifftshift(xGrid)/xsize;
function [yGrid, xGrid, zGrid] = freqgrid3(ysize,xsize, zsize)

%  [yGrid, xGrid, Grid] = freqgrid(ysize,xsizezsize)
%
% Creates a 3d grid of frequencies of a given size, respecting the MATLAB
% convention for the location of the DC component (top left). This means
% that if a filter is desgined on this grid, ifftn(filter) will behave as
% expected. Use fftshift before visualising.
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

% Dimensions
ymid = floor(ysize/2); 
xmid = floor(xsize/2);
zmid = floor(zsize/2); 

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

if(mod(zsize,2)) == 0
    zmax = zmid-1;
else
    zmax = zmid;
end

% Create the grid for the filter
[yGrid, xGrid, zGrid] = ndgrid( -ymid : ymax, -xmid : xmax,  -zmid : zmax);
yGrid = ifftshift(yGrid);
xGrid = ifftshift(xGrid);
zGrid = ifftshift(zGrid);

yGrid = yGrid ./ ysize;
xGrid = xGrid ./ xsize;
zGrid = zGrid ./ zsize;
function LO = localOrientation(m2,m3,degrees)
%  function LO = localOrientation(m2,m3,degrees = false)
%
% This function calculates the local image orientation from the odd parts
% of the monogenic signal. The value returned is in radians unless the
% third argument is set to true. The angle is measured anti-clockwise from
% the increasing x-axis and corresponds to the direction of increasing
% intensity in the image. The value is wrapped to be between -pi and pi
% radians (or -180 to 180 degrees).
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

if (nargin < 3)
    degrees = false;
end

if (degrees)
    LO = atan2d(-m2,m3);
else
    LO = atan2(-m2,m3);
end

function LP = localPhase(m1,m2,m3, wl_ind)
%
%    LP = localPhase(m1,m2,m3,wl_ind)
%
% Calculates the local phase from a monogenic signal (m1,m2,m3) of an image
% or video. The phase is calculated for each scale independently.
%
% Alternatively, wl_ind, a parameter selecting the index of the wavelength
% to use may be passed.
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

if nargin < 4
   LP = atan2( sqrt(m2.^2 + m3.^2), m1);
else
   LP = atan2( sqrt(m2(:,:,:,wl_ind).^2 + m3(:,:,:,wl_ind).^2), m1(:,:,:,wl_ind));
end





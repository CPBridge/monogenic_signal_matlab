function PC = phaseCongruency(m1,m2,m3,T)
%      PC = phaseCongruency(m1,m2,m3,T=0.0)
%
% Calculate the phase congruency measure from the monogenic signal
% (m1,m2,m3,m4). This algorithm is from Felsberg and Sommer "A New
% Extension of Linear Signal Processing for Estimating Local Properties and
% Detecting Features" and requires that the monogenic signal is calculated
% at exactly 2 scales.
%
% T is a constant noise threshold used to supress the phase congruency in
% parts of the image with low energy at both scales. It is expressed as a
% value in the range 0-1 where 0 will supress no noise and 1 will suppress
% the entire image.
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

if (nargin < 4)
    T = 0.0;
elseif( T < 0.0 || T > 1.0)
    error('T should lie in the range 0-1')
end

if(any([size(m1,4), size(m2,4), size(m3,4)] ~= 2 ))
   error('You must provide a monogenic signal at exactly two wavelengths')
end

% Stack the parts of the vectors together for scales 1 and 2
f1 = cat(4,m1(:,:,:,1),m2(:,:,:,1),m3(:,:,:,1));
f2 = cat(4,m1(:,:,:,2),m2(:,:,:,2),m3(:,:,:,2));

% Dot product of the two vectors
dotprod = sum(f1.*f2,4);

% Norm of the two vectors
norm1 = sqrt(sum(f1.^2,4));
norm2 = sqrt(sum(f2.^2,4));

% Cross product of the two vectors amd its norm
xprod = cross(f1,f2,4);
normxprod = sqrt(sum(xprod.^2,4));

% Phase congruency
PC = max(dotprod-T*max(dotprod(:)),0)./(normxprod+norm1.*norm2);

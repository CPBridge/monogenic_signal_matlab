% This is an example script for using the monogenic signal code for 2D
% images

% Add monogenic_signal source directory to path
addpath('src')

% Load a matlab test image
% If running under octave, you will need to change this to an image file of
% your choice.
% Note that the monogenic signal is intended to on greyscale images (using it on
% a colour image will result in the three channels being processed independently).
I = imread('coins.png');
[Y,X] = size(I);

% First we have to choose a set of centre-wavelengths for our filters,
% typically you will want to play around with this a lot.
% Centre-wavelengths are expressed in pixel units. Here we use a set of
% wavelenths with a constant scaling factor of 1.5 between them, starting
% at 20 pixels
cw = 20*1.5.^(0:4);

% Now use these wavelengths to create a structure containing
% frequency-domain filters to calculate the mnonogenic signal. We can
% re-use this structure many times if we need for many images of the same
% size and using the same wavelength. We can choose from a number of
% different filter types, with log-Gabor ('lg') being the default. For lg
% filters we can also choose the shape parameter (between 0 and 1), which
% governs the bandwidth (0.41 gives a three-octave filter, 0.55 gives a two
% octave filter)
filtStruct = createMonogenicFilters(Y,X,cw,'lg',0.55);

% Now we can use this structure to find the monogenic signal for the image
[m1,m2,m3] = monogenicSignal(I,filtStruct);

% The returned values are the three parts of the monogenic signal: m1 is
% the even part, and m2 and m3 are the odd parts in the vertical and
% horizontal directions respectively. Each array is Y x X x 1 x W, where
% X and Y are the image dimensions and W is the number of wavelengths.
% The filter responses to the filters of each scale are stacked along the
% fourth dimension.

% (Alternatively one may pass a 3D volume to monogenicSignal, in which case
% the 2D monogenic signal is found for each of the Z 2D planes independently and
% returned as a set of Y x X x Z x W arrays)

% From here we can straightforwardly find many of the derived measures by
% passing these three arrays

% Local energy (calculated on a per-scale basis)
LE = localEnergy(m1,m2,m3);

% Local phase (calculated on a per-scale basis)
LP = localPhase(m1,m2,m3);

% Local orientation (calculated on a per-scale basis)
% Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
LO = localOrientation(m2,m3);

% Feature symmetry and asymmetry (see Kovesi "Symmetry and Asymmetry from
% Local Phase") pick out 'blob-like' and 'boundary-like' structures
% respectively. This combines all the scales to give a single 2D image.
[FS,FA] = featureSymmetry(m1,m2,m3);

% Oriented symmetry and asymmetry are like the above but contain more
% information. Symmetric blobs are differentiated into peaks and troughs
% by the sign of the signed feature symmetry (SFS) measure. Oriented
% feature asymmetry (OFA) is a vector quantity (returned below as separate
% components) that describes both the magnitude and orientation of the
% boundary
[OFA_y,OFA_x,SFS] = orientedSymmetry(m1,m2,m3);

% Display
figure()
imshow(I), axis image, axis off, colormap gray
title('Test Image')


figure()
imagesc(reshape(LE,Y,[])), axis image, axis off, colormap gray
title('Local Energy Over Scales')

figure()
imagesc(reshape(LP,Y,[])), axis image, axis off, colormap gray
title('Local Phase Over Scales')

figure()
imagesc(reshape(LO,Y,[])), axis image, axis off, colorbar
title('Local Orientation Over Scales (radians)')

figure()
imagesc([FS,FA]), axis image, axis off, colormap gray
title('Feature Symmetry and Asymmetry')

figure()
imagesc([FS,OFA_y,OFA_x]), axis image, axis off, colormap gray
title('Signed Feature Symmetry and Oriented Feature Asymmetry')

%%

% If you want to visualise the filters to better understand what's going on,
% we can do something like this:
% First the frequency-domain representation (which is how the filter is stored
% and used) for just the third filter in the stack bpFilt(:,:,1,3)
figure()
subplot(1,3,1)
surf(real(fftshift(filtStruct.bpFilt(:,:,1,3))),'edgecolor','none')
title('Even Part of The Frequency Domain Filter')
subplot(1,3,2)
surf(real(fftshift(filtStruct.bpFilt(:,:,1,3).*filtStruct.ReiszFilt)),'edgecolor','none')
title('First Odd Part of The Frequency Domain Filter')
subplot(1,3,3)
surf(imag(fftshift(filtStruct.bpFilt(:,:,1,3).*filtStruct.ReiszFilt)),'edgecolor','none')
title('Second Odd Part of The Frequency Domain Filter')

% Now the image (spatial) domain filters, which can be found via ifft
figure()
subplot(1,3,1)
surf(real(fftshift(ifft2(filtStruct.bpFilt(:,:,1,3)))),'edgecolor','none')
title('Even Part of The Image Domain Filter')
subplot(1,3,2)
surf(real(fftshift(ifft2(filtStruct.bpFilt(:,:,1,3).*filtStruct.ReiszFilt))),'edgecolor','none')
title('First Odd Part of The Image Domain Filter')
subplot(1,3,3)
surf(imag(fftshift(ifft2(filtStruct.bpFilt(:,:,1,3).*filtStruct.ReiszFilt))),'edgecolor','none')
title('Second Odd Part of The Image Domain Filter')



%% Now for phase congruency
clear

% Load a matlab test image and convert to greyscale
I = rgb2gray(imread('board.tif'));
[Y,X] = size(I);

% This time we have to use exactly two scales, as this is required by
% Felsberg's phase congruency method. As he suggests, let's use
% three-ocatave filters and leave a three-octave spacing between them.
% We want small filters to pick out the fine detail in this image.
cw = [3,24];

% Construct new filters, as before
filtStruct = createMonogenicFilters(Y,X,cw,'lg',0.41);

% Find monogenic signal, as before
[m1,m2,m3] = monogenicSignal(I,filtStruct);

% Now use the phase congruency algorithm. The fourth parameter is a
% threshold between 0 and 1 used for noise supression. You will always need
% to use this to get reasonable results. Somewhere between 0 and 0.1 should
% do in most cases.
PC = phaseCongruency(m1,m2,m3,0.05);

% Display
figure()
imshow(I), axis image, axis off, colormap gray
title('Test Image')

figure()
imagesc(PC), axis image, axis off, colormap gray
title('Phase Congruency')

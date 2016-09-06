% This is an example script for using the monogenic signal code for 3D
% volume images. Most things are the same, with comparable function names
% with '3D' added

% Add monogenic_signal source directory to path
addpath('src')

% Load a matlab test volume image, the 'mri' dataset
% If running under Octave, you will need to provide a different volume image
% here
load mri
D = squeeze(D); % get rid of the third singleton dimension
[Y,X,Z] = size(D);

% First we have to choose a set of centre-wavelengths for our filters,
% For the 3D case you can choose a different set of wavelengths for the Z
% dimension from those used for the X-Y plane. This is because I was
% working with data that were not sampled isotropically. Normally you
% should stick to the same set of wavelengths (in fact the mri dataset is
% not sampled isotropically but we'll pretend it is...)
cw = 5*1.5.^(0:4);

% Now use these wavelengths to create a structure containing
% frequency-domain filters to calculate the monogenic signal. Same as for
% the 2D case, except now we need to pass three image dimensions and the
% two sets of wavelengths
filtStruct = createMonogenicFilters3D(Y,X,Z,cw,cw,'lg',0.55);

% Now we can use this structure to find the monogenic signal for the volume
[m1,m2,m3,m4] = monogenicSignal3D(D,filtStruct);

% Now we have four parts of the monogenic signal as there is another odd
% part for the Z direction. Each is an Y x X x Z x W array where again W is
% the number of scales/wavelengths used, i.e. the responses to each filter
% are again stacked along the fourth dimension

% From here we can straightforwardly find many of the derived measures by
% passing these four arrays

% Local energy (calculated on a per-scale basis)
LE = localEnergy3D(m1,m2,m3,m4);

% Local phase (calculated on a per-scale basis)
LP = localPhase3D(m1,m2,m3,m4);

% Note that we do not provide a local orientation function for 3D. This is
% because there are various conventions for describing orientation in 3D
% space and people's requirements will vary. It is straightforward to write
% your own for your requirements.

% Feature symmetry and asymmetry (see Kovesi "Symmetry and Asymmetry from
% Local Phase") pick out 'blob-like' and 'boundary-like' structures
% respectively. Combines all the scales to give a single 3D image.
[FS,FA] = featureSymmetry3D(m1,m2,m3,m4);

% Display one slice
slice = 13; % near the middle
figure()
imshow(D(:,:,slice)), axis image, axis off, colormap gray
title('Test Volume Slice')


figure()
imagesc(reshape(LE(:,:,slice,:),Y,[])), axis image, axis off, colormap gray
title('Local Energy Over Scales')

figure()
imagesc(reshape(LP(:,:,slice,:),Y,[])), axis image, axis off, colormap gray
title('Local Phase Over Scales')

figure()
imagesc([FS(:,:,slice),FA(:,:,slice)]), axis image, axis off, colormap gray
title('3D Feature Symmetry and Asymmetry')

%% Now for phase congruency
% We'll use the same test volume here

% This time we have to use exactly two scales, as this is required by
% Felsberg's phase congruency method. As he suggests, let's use
% three-ocatave filters and leave a three-octave spacing between them.
cw = [6,48];

% Construct new filters, as before
filtStruct = createMonogenicFilters3D(Y,X,Z,cw,cw,'lg',0.55);

% Find monogenic signal, as before
[m1,m2,m3,m4] = monogenicSignal3D(D,filtStruct);

% Now use the 3D phase congruency algorithm. Just like the 2D case
PC = phaseCongruency3D(m1,m2,m3,m4,0.05);

% Display
slice = 13; % near the middle
figure()
imshow(D(:,:,slice)), axis image, axis off, colormap gray
title('Test Volume Slice')

figure()
imagesc(PC(:,:,slice)), axis image, axis off, colormap gray
title('3D Phase Congruency')

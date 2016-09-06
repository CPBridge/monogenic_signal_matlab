# Monogenic Signal MATLAB (and GNU Octave) Implementation

This repository contains MATLAB code to calculate the monogenic signal (Felsberg and Sommer) for 2D and 3D images, as well as many quantities that can be derived from the monogenic signal such as Feature Symmetry and Asymmetry, and Phase Congruency.

The monogenic signal is an alternative way of representing an image, which has a number of advantages for further processing. For an introduction to the monogenic signal and derived features with references to the relevant scientific literature, please see [this document](https://chrispbridge.files.wordpress.com/2016/05/monogenic2.pdf) (PDF).

### Capabilities

Functions are provided to calculate the following quantities for 2D and 3D images:

* Monogenic Signal
* Local Energy, Local Phase and Local Orientation to describe the local properties of image
* Feature Symmetry and Asymmetry, respond to symmetric 'blobs' and boundaries with robustness to variable contrast
* Oriented Feature Symmetry and Asymmetry, as above but also containing the polarity of the symmetry and the orientatation of the boundaries
* Phase Congruency, responds to edges in the image

### Instructions For Use

Read through and run the heavily-commented example files (example2D.m and example3D.m) in order to learn how to use the functions (you will need to change the image files used in the examples if you are using Octave). Each function also has a relatively complete description that can be read in the accessed via the help interface

### Compatibility and Dependencies
MATLAB or GNU Octave (all versions on all operating systems should work as far as I am aware).

### Installation
You should not need to do anything to install except ensure that the directory containing the source files are in your MATLAB/Octave path whenever you want to use the it.

### Contributors
This software was written by Christopher Bridge (University of Oxford) and was based on previous code by Ana Namburete and Vicente Grau.

### Licence
This software is licensed under the GNU Public Licence (GPL). You are free to edit and distribute this code providing certain conditions are met. Read the full licence for further information.

### Publications

This code has been used in the following publication:

* C.P. Bridge and J.A. Noble, “[Object Localisation in Fetal Ultrasound Images Using Invariant Features](http://ieeexplore.ieee.org/document/7163839/)”. Proceedings of IEEE International Symposium on Biomedical Imaging, New York City, 2015

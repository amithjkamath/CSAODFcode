## CSAODFcode

This is a MATLAB implementation of the CSA ODF Algorithm for white matter fiber estimation using Diffusion MR data. More details in my [masters' thesis](http://conservancy.umn.edu/handle/11299/140183).

To demonstrate how this works, run the script:

`runExperiment.mlx`

This sets up the experiment parameters, creates dummy data to analyze, computes the Spherical Harmonics coefficients based on the gradient table, and computes the transformation matrix, and then processes the data.

## CLASSES

* createDiffusionData
* csaOdfExperiment
* gradientTable
* sphericalHarmonicsMatrix
* volume

## ACKNOWLEDGEMENTS

Depends on MATLAB and the Curve Fitting Toolbox from MathWorks Inc.
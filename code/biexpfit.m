function [fitresult, gof] = biexpfit(bvals,EE)
%    BIEXPFIT: this does a bi-exponential fit on the diffusion data. 
%    bvals is the array of bvalues on which the data is fitted on.
%    EE is the discrete measured values (measured or evaluated from the shell itself).
%    This function depends on the Curve Fitting Toolbox from MathWorks.

ft = fittype( 'exp2' );

opts = fitoptions( ft );
opts.Display = 'Off';
opts.Robust = 'LAR';

opts.Upper = [1 0 1 0];
opts.Lower = [0 -1 0 -1];

opts.MaxIter = 800; %Higher than the default 400.
opts.MaxFunEvals = 800;
%opts.StartPoint = [0.75 -0.005 0.25 -0.005];
%opts.Algorithm = 'Levenberg-Marquardt';

% Fit model to data.
[fitresult, gof] = fit( bvals, EE, ft, opts );

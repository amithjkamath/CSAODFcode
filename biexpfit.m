function [fitresult, gof] = biexpfit(bvals,EE)
ft = fittype( 'exp2' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 -1 0 -1];
opts.Robust = 'LAR';
opts.Upper = [1 0 1 0];
opts.MaxIter = 500; %Higher than the default 400.
%opts.StartPoint = [0.75 -0.005 0.25 -0.005];
%opts.Algorithm = 'Levenberg-Marquardt';

% Fit model to data.
[fitresult, gof] = fit( bvals, EE, ft, opts );

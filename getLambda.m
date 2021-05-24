function [lambda] = getLambda(fibers, bvalue, SNR, order)
    % lambda is a function of fibers, bvalue, SNR and order of reconstruction.

    % This assumes that bvalues are specified in order of 10^6 mm2/s.
    % The numbers calculated here are based on SNR divided by 5, and so this
    % conversion.

    % Author: Amith Kamath
    if nargin ~= 4
       error('Not enough input arguments. Fibers, bvalue, SNR and order are required.'); 
    end
    parseInputs(fibers, bvalue, SNR, order);

    % calculate the lambda value after these conversions. 
    
    bvalue = bvalue/1000;
    SNR = SNR/5;

    switch order
        case 4
            switch fibers
                case 1
                    lambda = [7.629 -0.9839 -2.694 0.04133 0.145 0.2613]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 2
                    lambda = [12.23 -1.603 -4.381 0.06918 0.2333 0.413]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 3
                    lambda = [16.16 -2.012 -5.832 0.08231 0.3006 0.5437]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
            end
        case 6
            switch fibers
                case 1
                    lambda = [11.16 -1.118 -4.149 0.04554 0.1671 0.3908]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 2
                    lambda = [15.5 -1.565 -5.765 0.06596 0.2307 0.5424]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 3
                    lambda = [19.13 -1.985 -7.012 0.08712 0.285 0.6544]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
            end
        case 8
             switch fibers
                case 1
                    lambda = [11.91 -1.113 -4.308 0.0469 0.1633 0.3956]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 2
                    lambda = [16.72 -1.684 -6.051 0.07119 0.2479 0.5538]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
                case 3
                    lambda = [19.29 -1.866 -6.926 0.07831 0.2736 0.634]*[1 SNR bvalue SNR.^2 SNR*bvalue bvalue.^2]';
             end
    end

    % lamda is scaled down due to the computations.
    lambda = lambda/1000;
    
    if(lambda < 0)
        error('Something went wrong. Lambda is negative. Please review the inputs to getLambda.');
    end
end

function parseInputs(fibers, bvalue, SNR, order)
    if(fibers ~= 1 && fibers ~= 2 && fibers ~= 3)
        error('Number of fibers must be 1, 2 or 3. Lambda is defined only for these values and no other.');
    end

    if length(bvalue) ~= 1
       error('bvalue must be a scalar value'); 
    end
    
    if(SNR < 3 || SNR > 50)
        error('SNR has to be between 3 and 50 for accurate results.');
    end
    
    if(bvalue < 500 || bvalue > 10000)
        error('bvalue must be between 500 and 10000. Other values are not tested for, and will be inaccurate.');
    end

    if(order ~= 4 && order ~= 6 && order ~= 8)
        error('SH reconstruction order must be 4, 6, or 8. Lambda is defined only for these values and no other.');
    end
end
classdef createDiffusionData < matlab.mixin.SetGet
    
    properties
        fibers = 1;
        bvals = [];
        SNR = 50;
        shellN = [];
        bvecs = [];
        data = [];
    end
    
    methods
        function diff = createDiffusionData()
            diff.fibers = 1;
            diff.bvals = [];
            diff.SNR = 50;
            diff.shellN = [];
            diff.bvecs = [];
        end
        
        function make(diff,fibersIn, bvalsIn, SNRIn, shellNIn, bvecsIn, anglesIn)
            S0_ = 250;
            nshells = max(shellNIn);
            nbvecs = size(bvecsIn, 1);
            E  = zeros(nbvecs, 1);
            S0 = diff.DTricedist(S0_, S0_/SNRIn, nbvecs);
            
            if (size(bvalsIn,1) == nshells)
                Q1 = qr(eye(3));
                %Q1 = qr(rand(3));
                Q1 = Q1./det(Q1);
                
                %Fix the angle of the fibers here. Choose Rand to simulate random fibers
                %upto rnge specified by the multiplier.
                
                angle = anglesIn;
                
                Amp = [1.7 0 0; 0 0.2 0; 0 0 0.2];
                
                D1 = Q1*Amp*Q1'.* 10^-3;
                
                axis1 = Q1(:,1)';
                Q2 = diff.rotMatrix(axis1,(pi/2) - angle)*Q1;
                
                D2 = Q2*[0.2 0 0; 0 1.7 0; 0 0 0.2]*Q2'.* 10^-3;
                
                axis2 = Q2(:,2)';
                Q3 = diff.rotMatrix(axis2,(pi/2)+angle)*Q2;
                
                D3 = Q3*[0.2 0 0; 0 0.2 0; 0 0 1.7]*Q3'.* 10^-3;
                
                % Create DWI signal
                for i=1:nbvecs
                    u = squeeze(bvecsIn(i,:))';
                    bvalue = bvalsIn(shellNIn(i)); % This is where the bvalue is selected from the shell number specified in the data.
                    switch fibersIn
                        case 1
                            E(i) = S0(i) * (exp(-bvalue*u'*D1*u)); % j indexes the number of bvalues for that particular direction, along i.
                        case 2
                            E(i) = S0(i) * (0.6*exp(-bvalue*u'*D1*u) + 0.4*exp(-bvalue*u'*D3*u)); % j indexes the number of bvalues for that particular direction, along i.
                        case 3
                            E(i) = S0(i) * (0.4*exp(-bvalue*u'*D1*u) + 0.3*exp(-bvalue*u'*D3*u) + 0.3*exp(-bvalue*u'*D2*u)); % j indexes the number of bvalues for that particular direction, along i.
                    end
                end
                data3X = E;
                B0 = diff.DTricedist(S0_, S0_/SNRIn, 1);
                diff.data = data3X./B0;
            end
        end
    end
    
    methods (Access = private)
        function R = DTricedist(~, devia, sigma, numpoints)
            theta = rand(1);
            x = (sigma).*randn(numpoints,1) + devia*cos(theta);
            y = (sigma).*randn(numpoints,1) + devia*sin(theta);
            R = sqrt(x.^2 + y.^2);
        end
        
        function outM = rotMatrix(~, axis, angle)
            % axis has to be a row vector for this to work.
            % angle in radians please!
            % Rotation is clockwise.
            outM = eye(3).*cos(angle) + sin(angle).*[0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0] + (1 - cos(angle)).*(axis'*axis);
            % outM = (1/sqrt(2)).*outM;
        end
        
    end
end
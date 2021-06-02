classdef csaodfExperiment < matlab.mixin.SetGet
    
    properties
        fibers = 1;
        bvals = [];
        SNR = 40;
        order = 8;
        lambda = [];
    end
    
    methods
        function expt = csaodfExperiment()
            expt.fibers = 1;
            expt.bvals = [];
            expt.SNR = 5;
            expt.order = 4;
            expt.lambda = [];
        end
        
        function set.fibers(expt, inFibers)
            if(inFibers ~= 1 && inFibers ~= 2 && inFibers ~= 3)
                error('Number of fibers must be either 1, 2 or 3.');
            else
                expt.fibers = inFibers;
            end
        end
        
        function set.order(expt, inOrder)
            if(inOrder ~= 4 && inOrder ~= 6 && inOrder ~= 8)
                error('Order has to be either 4, 6 or 8. No other values are allowed!');
            else
                expt.order = inOrder;
            end
        end
        
        function set.SNR(expt, inSNR)
            if(inSNR < 5 || inSNR > 50)
                error('SNR has to be between 5 and 50. No other values are allowed!');
            else
                expt.SNR = inSNR;
            end
        end
        
        function vol = processVolume(expt, v, gt, sh)
            
            if(~isa(v,'volume'))
                error('First input value must be of class volume.');
            end
            if(~isa(gt,'gradientTable'))
                error('Second input value must be of class gradientTable.');
            end
            if(~isa(sh,'sphericalHarmonicsMatrix'))
                error('Third input value must be of class sphericalHarmonicsMatrix.');
            end
            
            numSH = (expt.order + 1)*(expt.order + 2)/2;
            vol = volume;
            vol.imageData = zeros(v.getxdim, v.getydim, v.getzdim, numSH);
            
            for i = 1:v.getxdim
                for j = 1:v.getydim
                    for k = 1:v.getzdim
                        % Do the data interpolation on the new bvecs directions.
                        % gt.table is passed in twice in the function below
                        % because the 'full' estimates are made in the same
                        % directions as earlier.
                        dataInt = expt.interpolateData(v.dataAtVoxel(i,j,k,0,0), gt.table, gt.table, gt.shellInd, expt.order, expt.lambda);
                        
                        %Append b0 values to the data. This is done after the interpolation on the new directions, which does not hence interfere with the reconstruction accuracy.
                        b0dataInt = [ones(size(gt.table,1),1), dataInt];
                        
                        fprintf('Processing voxel (  %d,  %d,  %d) ..',i,j,k);
                        %% This is the biexponential fit using the Trust region method, Matlab.
                        %Preallocate!
                        dataC1 = zeros(size(gt.table,1),1);
                        dataC2 = zeros(size(gt.table,1),1);
                        dataex1 = zeros(size(gt.table,1),1);
                        dataex2 = zeros(size(gt.table,1),1);
                        dataMod = zeros(size(gt.table,1),1);
                        
                        %% Do the biexponential fit.
                        for x = 1:size(gt.table,1)
                            if(mod(x,size(gt.table)/20)==0)
                                fprintf('.');
                            end
                            fitresult = biexpfit([0; gt.bValues], b0dataInt(x,:)');
                            dataC1(x) = fitresult.a;
                            dataC2(x) = fitresult.c;
                            
                            dataex1(x) = fitresult.b;
                            dataex2(x) = fitresult.d;
                            dataMod(x) = dataC1(x)*log(-(dataex1(x))) + dataC2(x)*log(-dataex2(x));
                        end
                        
                        % This is a failure mode - don't use lambda in the
                        % case where the optimization hasn't converged.
                        bvecsN = gt.table(dataC1 ~= 0.5,:);
                        dataModN = dataMod(dataC1 ~= 0.5);
                        
                        if(numel(bvecsN) ~= numel(gt.table))
                            shN = sphericalHarmonicsMatrix(bvecsN, expt.order, 0);
                            SH = shN.computeODF(dataModN);
                        else 
                            SH = sh.computeODF(dataModN);
                        end
                        
                        vol.imageData(i,j,k,:) = SH;
                        fprintf(' Done!\n');
                    end
                end
            end
        end
        
        function visualizeODF(~, dataPoints, directions)
            % to do.
            fcs = convhulln(directions);
            
            points(:,1) = dataPoints.*directions(:,1);
            points(:,2) = dataPoints.*directions(:,2);
            points(:,3) = dataPoints.*directions(:,3);
            
            patch('Vertices',[points(:,1) points(:,2) points(:,3)],'Faces',fcs,'FaceVertexCData',dataPoints,...
                'FaceColor','interp','EdgeColor','interp','BackFaceLighting','lit',...
                'FaceLighting','phong');
            grid on;
            xlabel('X axis');
            ylabel('Y axis');
            zlabel('Z axis');
            view([1 1 1]);
            colorbar;
        end
    
        function points = genNSpherePoints(~, n, varargin)
            if nargin == 2
                r = 1;
            else
                assert(nargin <= 3)
                r = varargin{2};
            end
            assert(n > 0)
            theta = 2.*pi.*rand([n 1]);
            phi   = acos(1 - 2.*rand([n 1]));
            x     = r.*sin(phi).*cos(theta);
            y     = r.*sin(phi).*sin(theta);
            z     = r.*cos(phi);
            points = [x y z];
        end
        
        function computeLambda(expt)
            % lambda is a function of fibers, bvalue, SNR and order of reconstruction.
        
            % This assumes that bvalues are specified in order of 10^6 mm2/s.
            % The numbers calculated here are based on SNR divided by 5, and so this
            % conversion.
        
            % Author: Amith Kamath
            if(expt.fibers ~= 1 && expt.fibers ~= 2 && expt.fibers ~= 3)
                error('Number of fibers must be 1, 2 or 3. Lambda is defined only for these values and no other.');
            end
            
            if(expt.SNR < 3 || expt.SNR > 50)
                error('SNR has to be between 3 and 50 for accurate results.');
            end
        
            if(expt.order ~= 4 && expt.order ~= 6 && expt.order ~= 8)
                error('SH reconstruction order must be 4, 6, or 8. Lambda is defined only for these values and no other.');
            end
            for idx = 1:numel(expt.bvals)
                if length(expt.bvals(idx)) ~= 1
                   error('bvalue must be a scalar value'); 
                end
                
                if(expt.bvals(idx) < 1000 || expt.bvals(idx) > 5000)
                    error('bvalue must be between 1000 and 5000. Other values are not tested for, and will be inaccurate.');
                end
                % calculate the lambda value after these conversions. 
                
                iBvalue = expt.bvals(idx)/1000;
                iSNR = expt.SNR/5;
            
                switch expt.order
                    case 4
                        switch expt.fibers
                            case 1
                                expt.lambda(idx) = [7.629 -0.9839 -2.694 0.04133 0.145 0.2613]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 2
                                expt.lambda(idx) = [12.23 -1.603 -4.381 0.06918 0.2333 0.413]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 3
                                expt.lambda(idx) = [16.16 -2.012 -5.832 0.08231 0.3006 0.5437]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                        end
                    case 6
                        switch expt.fibers
                            case 1
                                expt.lambda(idx) = [11.16 -1.118 -4.149 0.04554 0.1671 0.3908]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 2
                                expt.lambda(idx) = [15.5 -1.565 -5.765 0.06596 0.2307 0.5424]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 3
                                expt.lambda(idx) = [19.13 -1.985 -7.012 0.08712 0.285 0.6544]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                        end
                    case 8
                         switch expt.fibers
                            case 1
                                expt.lambda(idx) = [11.91 -1.113 -4.308 0.0469 0.1633 0.3956]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 2
                                expt.lambda(idx) = [16.72 -1.684 -6.051 0.07119 0.2479 0.5538]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                            case 3
                                expt.lambda(idx) = [19.29 -1.866 -6.926 0.07831 0.2736 0.634]*[1 iSNR iBvalue iSNR.^2 iSNR*iBvalue iBvalue.^2]';
                         end
                end
            
                % lamda is scaled down due to the computations.
                expt.lambda(idx) = expt.lambda(idx)/1000;
                
                if(expt.lambda(idx) < 0)
                    error('Something went wrong. Lambda is negative. Please review the inputs to getLambda.');
                end
            end
        end
    end
    
    methods (Access = private)
        function dataInt = interpolateData(~, data3X, bvecsI, bvecsO, shellN, order, lambda)
            nshells = max(shellN);
            
            for ns = 1:nshells
                bvecsSH = bvecsI((shellN == ns),:);
                dataSH = data3X(shellN == ns);
                sh = sphericalHarmonicsMatrix(bvecsSH, order, lambda(ns));
                sigSH(ns,:) = sh.sampleToCoeff(dataSH);
            end
            
            shO = sphericalHarmonicsMatrix(bvecsO, order, 0);
            
            for ns = 1:nshells
                dataInt(ns,:) = shO.T*sigSH(ns,:)';
            end
            
            dataInt = dataInt';
        end
    end
end


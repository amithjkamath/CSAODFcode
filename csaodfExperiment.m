classdef csaodfExperiment < matlab.mixin.SetGet
    % csaodfExperiment
    % Summary of this class goes here
    %   Detailed explanation goes here
    
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
            expt.SNR = 40;
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
        
        function set.bvals(expt, inBvals)
            % no input validation as of now.
            expt.bvals = inBvals;
        end
        
        function set.lambda(expt, inLambda)
            % no input validation as of now.
            expt.lambda = inLambda;
        end
        
        function vol = processVolume(expt, v, gt, sh)
            % to do.
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
                        dataInt = expt.interpolateData(v.getValues(i,j,k,0,0), gt.table, gt.table, gt.shellInd, expt.order, expt.lambda);
                        
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
                        bvecsN = gt.table(dataC1 ~= 0.5,:);
                        dataModN = dataMod(dataC1 ~= 0.5);
                        
                        if(bvecsN ~= gt.table)
                            sh.make(bvecsN,expt.order,0);
                        end
                        
                        SH = sh.getODF(dataModN);
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
    end
    
    methods (Access = private)
        function dataInt = interpolateData(~, data3X, bvecsI, bvecsO, shellN, order, lambda)
            nshells = max(shellN);
            
            for ns = 1:nshells
                bvecsSH = bvecsI((shellN == ns),:);
                dataSH = data3X(shellN == ns);
                sh = sphericalHarmonicsMatrix;
                sh.make(bvecsSH,order,lambda(ns));
                sigSH(ns,:) = sh.getCoeff(dataSH);
            end
            
            shO = sphericalHarmonicsMatrix;
            shO.make(bvecsO,order,0);
            
            for ns = 1:nshells
                dataInt(ns,:) = shO.T*sigSH(ns,:)';
            end
            
            dataInt = dataInt';
        end
    end
end


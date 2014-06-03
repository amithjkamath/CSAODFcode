%% clear the workspace.
clear all
clc

%% Generate the data with the specified settings.

fibers = 2;
bvals = [1000 2000 3000];
SNR = 45;
order = 8;
volSize = 2;

gt = gradientTable;
gt.readFromTxt('samples3SH.txt',bvals);

v = volume;
diff = createDiffusionData;
diff.make(fibers, gt.bValues, SNR, gt.shellInd, gt.table, pi/2);
v.imageData = zeros(volSize,volSize,volSize,300);
v.imageData(1,1,1,:) = diff.data;
v.imageData(1,1,2,:) = diff.data;
v.imageData(1,2,1,:) = diff.data;
v.imageData(1,2,2,:) = diff.data;
v.imageData(2,1,1,:) = diff.data;
v.imageData(2,1,2,:) = diff.data;
v.imageData(2,2,1,:) = diff.data;
v.imageData(2,2,2,:) = diff.data;

SHv = volume;
SHv.imageData = zeros(volSize,volSize,volSize,15);

for x = 1:size(gt.bValues,1)
    lambda(x) = getLambda(fibers, gt.bValues(x),SNR,order);
end

%% Process Volume
for i = 1:volSize
    for j = 1:volSize
        for k = 1:volSize
            % Do the data interpolation on the new bvecs directions.
            dataInt = IntData(v.getValues(i,j,k,0,0), gt.table, gt.table, gt.shellInd, order, lambda);
            
            %Append b0 values to the data. This is done after the interpolation on the new directions, which does not hence interfere with the reconstruction accuracy.
            b0dataInt = [ones(size(gt.table,1),1), dataInt];
            
            %% This is the biexponential fit using the Trust region method, Matlab.
            %Preallocate!
            dataC1 = zeros(size(gt.table,1),1);
            dataC2 = zeros(size(gt.table,1),1);
            dataex1 = zeros(size(gt.table,1),1);
            dataex2 = zeros(size(gt.table,1),1);
            dataMod = zeros(size(gt.table,1),1);
            
            %% Do the biexponential fit.
            for x = 1:size(gt.table,1)
                fprintf('Fitting for %d out of %d directions.\n',x,size(gt.table,1));
                [fitresult, gof] = biexpfit([0; gt.bValues],b0dataInt(x,:)');
                dataC1(x) = fitresult.a;
                dataC2(x) = fitresult.c;
                
                dataex1(x) = fitresult.b;
                dataex2(x) = fitresult.d;
                dataMod(x) = dataC1(x)*log(-(dataex1(x))) + dataC2(x)*log(-dataex2(x));
            end
            
            %% Visualize the effective Diffusion profile using the biexp. model.
            % Optional: Remove the directions where the biexp. model has not converged.
            % This will make for a better looking output.
            bvecsN = gt.table(dataC1 ~= 0.5,:);
            dataModN = dataMod(dataC1 ~= 0.5);
            % figure, visODF(dataModN',bvecsN);
            
            %% display the interpolated data.
            %figure, visODF(dataModN, bvecsN);
            
            %% Create the ODF
            sh = SHmat;
            sh.make(bvecsN,4,0);
            SH = sh.getODF(dataModN);
            SHv.imageData(i,j,k,:) = SH;
        end
    end
end
%% Display the ODF
shDisp = SHmat;
shDisp.make(gt.table,4,0);
dataVals = shDisp.getData(SH);
figure, visODF(dataVals, gt.table);
view([0 1 0]);


fprintf('Setting up the experiment... ');
s = csaodfExperiment;

s.fibers = 2;
s.bvals = [1000, 2000, 3000];
s.SNR = 45;
s.order = 4;
% order 8 is way too descriptive. 4 is good!
for x = 1:size(s.bvals,2)
    s.lambda(x) = getLambda(s.fibers,s.bvals(x),s.SNR,s.order);
end

gt = gradientTable;
gt.readFromTxt('samples3SH.txt',s.bvals);
fprintf('Done! \n');

fprintf('Creating dummy data... ');
v = volume;
diffData = createDiffusionData;
diffData.make(s.fibers, gt.bValues, s.SNR, gt.shellInd, gt.table, pi/2);
v.imageData = zeros(1,1,2,300);
v.imageData(1,1,1,:) = diffData.data;
v.imageData(1,1,2,:) = diffData.data;
fprintf('Done! \n');

fprintf('Computing transformation matrix... ');
SH = sphericalHarmonicsMatrix;
SH.make(gt.table, s.order, 0);
fprintf('Done! \n');

fprintf('Processing data... \n');
ODFdata = s.processVolume(v, gt, SH);
fprintf('Done! \n');

fprintf('Creating visualization... ');
dataVals = SH.getData(ODFdata.getValues(1,1,1,0,0));
s.visualizeODF(dataVals,gt.table);
fprintf('Done! \n');
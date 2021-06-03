classdef testGradientTable < matlab.unittest.TestCase
    
    properties
        
    end
    
    methods(Test)
        
        function constructor(testCase)
            
            % Good inputs
            directions = rand(300, 3);
            shellIdx = randi(3, [300 1]);
            bVals = [1000 2000 3000];
            gradientTable(directions, shellIdx, bVals);
            
            % bad inputs
            testCase.verifyError(@()gradientTable([], shellIdx, bVals), "gradientTable:directionsMissing");
            
            testCase.verifyError(@()gradientTable(), "gradientTable:numArgsWrong");
            
            testCase.verifyError(@()gradientTable(directions), "gradientTable:numArgsWrong");
            
            testCase.verifyError(@()gradientTable(directions, shellIdx), "gradientTable:numArgsWrong");
            
            directions = rand(299, 3);
            testCase.verifyError(@()gradientTable(directions, shellIdx, bVals), "gradientTable:directionsMismatchShells");
            
            shellIdx = randi(3, [299 1]);
            bVals = [1000 2000];
            testCase.verifyError(@()gradientTable(directions, shellIdx, bVals), "gradientTable:bvalsMismatchShells");
        end
        
        function unnormalizedTable(testCase)
            directions = rand(300, 3);
            directions(1, :) = [1.5, 0.99, 0.98];
            shellIdx = randi(3, [300 1]);
            bVals = [1000 2000 3000];
            testCase.verifyError(@()gradientTable(directions, shellIdx, bVals), "gradientTable:unnormalizedTable");
        end
        
        function negativeShellIndex(testCase)
            directions = rand(300, 3);
            shellIdx = randi(2, [300 1]);
            shellIdx(255) = 0;
            bVals = [1000 2000 3000];
            testCase.verifyError(@()gradientTable(directions, shellIdx, bVals), "gradientTable:negativeShellIndex");
        end
        
        function negativeBValue(testCase)
            directions = rand(300, 3);
            shellIdx = randi(3, [300 1]);
            bVals = [1000 2000 -3000];
            testCase.verifyError(@()gradientTable(directions, shellIdx, bVals), "gradientTable:negativeBValue");            
        end
        
        function readFromTxtFile(testCase)
            % Positive test
            bvals = [1000, 2000, 3000];
            gradientTable.readFromTxt("samples3SH.txt", bvals);
            
            % Negative test
            testCase.verifyError(@()gradientTable.readFromTxt("emptyfile.txt", bvals), "gradientTable:txtNoDataInFile");                        
        end
        
        function readFromBvecBval(testCase)
            % Positive test
            gradientTable.readFromBvecBval("data/bvec.txt", "data/bvalue.txt");
            
            % Negative test
            testCase.verifyError(@()gradientTable.readFromBvecBval("emptyfile.txt", "emptyfile.txt"), "gradientTable:bvecbvalNoDataInFile");                         
        end
        
        function show(testCase)
            % Positive test
            directions = rand(300, 3);
            shellIdx = randi(3, [300 1]);
            bVals = [1000 2000 3000];
            gt = gradientTable(directions, shellIdx, bVals); 
            
            gt.show()
            close all
        end
    end
end
classdef gradientTable < matlab.mixin.SetGet
    %GRADIENTTABLE: class to hold data and allow transformation on the gradient table for Diffusion Tensor Imaging acquisition.
    %   Constructor:
    %   gt = gradientTable;
    %   empty gradient table. You can then initialize data from many different types of sources: scheme files, bval-bvec files or just text files.
    %   Properties include:
    %   table: 2D matrix including the x, y, z coordinates of the points where diffusion data is acquired.
    %   shellInd: If data is acquired on multiple Q-shells, this indicator array includes the index of shells.
    %   bValues: again, for multiple shell acquisitions, what bValues are they acquired at, is stored here. 
    %   for example, if acquisition is on three shells at 1000, 2000, 3000 units, bValues is [1000 2000 3000]. 
    %   This, along with shellInd provides a complete picture of what sample is acquired at what shell.
    %   Methods include:
    %   getter for bVecs (table) and read/write functions to scheme files, bval-bvec files and text files. 
    properties
        table = [];
        shellInd = [];
        bValues = [];
    end
    
    methods
        
        function GT = gradientTable(directions, shellIdx, bVals)
            if nargin ~= 3
                error("gradientTable:numArgsWrong", "Specify all arguments for gradientTable constructor")
            end
            if ~isempty(directions)
                % directions should be M-by-3
                assert(size(directions, 2) == 3)
                GT.table = directions;
                if( numel(bVals) ~= numel(unique(shellIdx)))
                    error("gradientTable:bvalsMismatchShells", "Number of bvalues must match number of shells")
                end
                if( size(directions, 1) ~= size(shellIdx, 1))
                    error("gradientTable:directionsMismatchShells", "Number of directions does not match shell indices")
                end
                GT.shellInd = shellIdx(:);
                GT.bValues = bVals(:);
            else
                error("gradientTable:directionsMissing", "First input should be M-by-3 and non-empty")
            end
        end
        
        function set.table(GT,tableIn)
            if((any(tableIn(:) > 1) || (any(tableIn(:) < -1))))
                error("gradientTable:unnormalizedTable", "Gradient table must be normalized. Incorrect numbers here.");
            end
            GT.table = tableIn;
        end
        
        function set.shellInd(GT,shellIdx)
            if(any(shellIdx(:) <= 0))
                error("gradientTable:negativeShellIndex", "Only positive indices allowed.");
            end
            GT.shellInd = shellIdx(:);
        end
        
        function set.bValues(GT,bVals)
            if(any(bVals(:) < 0))
                error("gradientTable:negativeBValue", "Only positive bvalues allowed.");
            end
            GT.bValues = bVals(:);
        end
        
        function show(GT)
            figure, hold on
            pointColor = ['r', 'g', 'b'];
            for i = 1:numel(GT.bValues)
                pointsInShell = GT.bValues(i).*GT.table(GT.shellInd == i, :);
                plot3(pointsInShell(:, 1), pointsInShell(:, 2), pointsInShell(:, 3), 'o', 'Color', pointColor(i))
            end
            hold off
            view(45, 45)
        end
    end
    
    methods (Static)

        function GT = readFromTxt(fileName, bVals)
            % txt files do not have bvalue info. They only have indices for
            % the shell that the sample lives in.
            dat = importdata(fileName);
            if ~isempty(dat)
                directions = dat.data(:,2:4);
                shellInd = dat.data(:,1);
                bValues = bVals(:);
                if( numel(bValues) ~= numel(unique(shellInd)))
                    error("gradientTable:txtBValsDontMatch", "Number of bvalues must match number of shells");
                end
                
                GT = gradientTable(directions, shellInd, bValues);
            else
                error("gradientTable:txtNoDataInFile", "Empty txt file. No data saved in gradientTable.");
            end
        end
        
        function GT = readFromBvecBval(bVecFile, bValFile)
            
            bvecs = importdata(bVecFile);
            bvals = importdata(bValFile);
            
            % Handle case of transposed entries.
            if size(bvals, 1) == 1
                bvals = bvals';
                bvecs = bvecs';
            end
            
            if ~isempty(bvecs) && ~isempty(bvals)
                if (size(bvecs,1) ~= size(bvals,1))
                    error("gradientTable:bvecbvalSampleMismatch", "Number of samples in bvecs does not match ones in bvals");
                end
                directions = bvecs(2:end, :);
                bValues = unique(bvals(2:end));
                shellInd = ones(size(directions,1), 1); 
                for i = 1:numel(bValues)
                    shellInd(bvals(2:end) == bValues(i)) = i;
                end
                if( numel(bValues) ~= numel(unique(shellInd)))
                    error("gradientTable:bvecbvalBvalShellMismatch", "Number of bvalues must match number of shells");
                end
                
                GT = gradientTable(directions, shellInd, bValues);
            else
                error("gradientTable:bvecbvalNoDataInFile", "Empty txt file. No data saved in gradientTable.");
            end
        end
    end
end


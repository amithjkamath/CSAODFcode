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
            if nargin == 0
                GT.table = [];
                GT.shellInd = [];
                GT.bValues = [];
                return;
            end
            if nargin ~= 3
                error('Specify all arguments for gradientTable constructor');
            end
            if ~isempty(directions)
                GT.table = directions;
                if( numel(bVals) ~= numel(unique(shellIdx)))
                    error('Number of bvalues must match number of shells');
                end
                GT.shellInd = shellIdx(:);
                GT.bValues = bVals(:);
            else
                GT.table = [];
                GT.shellInd = [];
                GT.bValues = [];
            end
        end
        
        function set.table(GT,tableIn)
            if((any(tableIn(:) > 1) || (any(tableIn(:) < -1))))
                error('Gradient table must be normalized. Incorrect numbers here.');
            end
            GT.table = tableIn;
        end
        
        function set.shellInd(GT,shellIdx)
            if(any(shellIdx(:) <= 0))
                error('Only positive indices allowed.');
            end
            GT.shellInd = shellIdx(:);
        end
        
        function set.bValues(GT,bVals)
            if(any(bVals(:) < 0))
                error('Only positive bValues allowed.');
            end
            GT.bValues = bVals(:);
        end
        
        function readFromScheme(GT,fileName)
            % scheme files have bvalue info as well.
            dat = importdata(fileName);
            if ~isempty(dat)
                GT.table = dat(:,1:3);
                GT.bValues = unique(dat(:,4));
                GT.shellInd = ones(size(dat,1),1);
                for i = 1:numel(GT.bValues)
                    GT.shellInd(dat(:,4) == GT.bValues(i)) = i;
                end
                if( numel(GT.bValues) ~= numel(unique(GT.shellInd)))
                    error('Number of bvalues must match number of shells');
                end
            else
                error('Empty scheme file. No data saved in gradientTable.');
            end
        end
        
        function readFromTxt(GT,fileName, bVals)
            % txt files do not have bvalue info. They only have indices for
            % the shell that the sample lives in.
            dat = importdata(fileName);
            if ~isempty(dat)
                GT.table = dat.data(:,2:4);
                GT.shellInd = dat.data(:,1);
                GT.bValues = bVals(:);
                if( numel(GT.bValues) ~= numel(unique(GT.shellInd)))
                    error('Number of bvalues must match number of shells');
                end
            else
                error('Empty txt file. No data saved in gradientTable.');
            end
        end
        
        function readFromBvecBval(GT,bVecFile, bValFile)
            bvecs = importdata(bVecFile)';
            bvals = importdata(bValFile)';
            if ~isempty(bvecs) && ~isempty(bvals)
                if (size(bvecs,1) ~= size(bvals,1))
                    error('Number of samples in bvecs does not match ones in bvals');
                end
                GT.table = bvecs;
                GT.bValues = unique(bvals(:));
                GT.shellInd = ones(size(bvals,1),1);
                for i = 1:numel(GT.bValues)
                    GT.shellInd(bvals == GT.bValues(i)) = i;
                end
                if( numel(GT.bValues) ~= numel(unique(GT.shellInd)))
                    error('Number of bvalues must match number of shells');
                end
            else
                error('Empty txt file. No data saved in gradientTable.');
            end            
        end
        
        function writeToSchemeFile(GT, outFileName)
            if ~isempty(GT.table)
                fid = fopen(outFileName, 'w');
                fprintf(fid, 'VERSION: 1\n');
                for i = 1:size(GT.table,1)
                    gamma = 2*pi*42.576*10^6;
                    bval = GT.bValues(GT.shellInd(i));

                    delta = 0.02;
                    G = 0.05;
                    TE = 0.1;
                    Delta = (bval/((gamma*delta*G)^2)) + (delta/3);
                    GT.table(i,:) = GT.table(i,:)./norm(GT.table(i,:));
                    %The *10^6 in the bvalues is to conform to the camino standard of using SI units.
                    fprintf(fid, '%+1.5f\t%+1.5f\t%+1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n', GT.table(i,1),GT.table(i,2),GT.table(i,3),G,Delta,delta,TE);
                end
                fclose(fid);
            else
                error('No data to be written to file.');
            end
        end
        
        function bVecs = getbvecs(GT,shellNum)
            if(shellNum < 0)
                error('shellNum needs to be positive.');
            elseif(shellNum == 0)
                bVecs = GT.table;
            else
                bVecs = GT.table(GT.shellInd == shellNum,:);
            end
        end
    end
    
end


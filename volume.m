classdef volume < hgsetget
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imageData = [];
        is4D = 0;
    end
    
    methods
        %setter for imageData
        function set.imageData(Vol,imageMat)
            if (ndims(imageMat) == 3)
                Vol.imageData = imageMat;          
            elseif (ndims(imageMat) == 4)
                Vol.imageData = imageMat;
                Vol.is4D = 1;
            elseif isempty(imageMat)
                Vol.imageData = [];                
            else
                error('Volume accepts only 3D and 4D data.');
            end
        end
        
        function set.is4D(Vol,value)
            if(value ~= 0 || value ~= 1)
                Vol.is4D = value;
            else
                error('is4D is either 0 or 1');
            end
        end
        
        % constructor.
        function Vol = volume(imageMat)
            if nargin == 1
                Vol.imageData = imageMat;
                if ndims(imageMat) == 4
                    Vol.is4D = 1;
                end
            elseif nargin == 0
                Vol.imageData = [];
            else
                error('Expecting one or no arguments.');
            end
        end
        
        function val = getxdim(Vol)
            val = size(Vol.imageData,1);
        end
        
        function val = getydim(Vol)
            val = size(Vol.imageData,2);
        end
        
        function val = getzdim(Vol)
            % this and getNumSamples has this logic because MATLAB returns
            % 1 even if there is no such dimension of data. zeros(100,100)
            % has a size(*,3) = 1, which I think is absurd.
            if size(Vol.imageData,3) == 1
                val = 0;
            else
                val = size(Vol.imageData,3);
            end
        end        
        
        function val = getNumSamples(Vol)
            if size(Vol.imageData,4) == 1
                val = 0;
            else
                val = size(Vol.imageData,4);
            end
        end
        
        %read from NIfTI file.
        function readFromNii(Vol, fileName)
            niiStruct = load_nii(fileName);
            Vol.imageData = niiStruct.img;
        end
        
        %write a NIfTI file.
        function writeToNii(Vol,fileName)
            error('No support for this for now. ');
        end
        
        function dataVals = getValues(Vol, inX, inY, inZ, shellIdx, shellNum)
            if isempty(Vol.imageData)
                error('No data in this volume. Cannot extract any!');
            end
            if(shellNum < 0)
                error('shellNum must be a positive number.');
            elseif (shellNum == 0 && shellIdx == 0)
                dataVals = squeeze(Vol.imageData(inX,inY,inZ,:));
            else
                dataVals = squeeze(Vol.imageData(inX,inY,inZ,shellIdx == shellNum));
            end
        end
    end
    
end


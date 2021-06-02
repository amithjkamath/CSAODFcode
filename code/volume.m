classdef volume < matlab.mixin.SetGet
    %VOLUME holds the 4 D volume of data. This is used to hold both the raw data and the processed SH coefficients.
    %   Usage notes:
    %   v = volume; %empty volume just constructed.
    %   public accessible imageData, accepts either 3D or 4D matrices only.
    %   properties include the raw 4D volume, and a boolean flag to check if it is 4D or not.
    %   methods include getting dimensions for calculation, reading and writing NIfTI files, 
    %   and a custom getter for accessing the element values in the 4th dimension.
    
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
        
        function val = getxdim(vol)
            val = size(vol.imageData, 1);
        end
        
        function val = getydim(vol)
            val = size(vol.imageData, 2);
        end
        
        function val = getzdim(vol)
            val = size(vol.imageData, 3);
        end        
        
        function val = getNumSamples(vol)
            val = size(vol.imageData, 4);
        end
        
        %read from NIfTI file.
        function readFromNii(Vol, fileName)
            imdata = niftiread(fileName);
            
            % First gradient direction is assumed to be B0. 
            Vol.imageData = imdata(:, :, :, 2:end)./imdata(:, :, :, 1);
            Vol.imageData = double(Vol.imageData);
        end
        
        function dataVals = dataAtVoxel(Vol, inX, inY, inZ, shellIdx, shellNum)
            if isempty(Vol.imageData)
                error('No data in this volume. Cannot extract any!');
            end
            if(shellNum < 0)
                error('shellNum must be a positive number.');
            elseif (shellNum == 0 && shellIdx == 0)
                dataVals = squeeze(Vol.imageData(inX, inY, inZ, :));
            else
                dataVals = squeeze(Vol.imageData(inX, inY, inZ, shellIdx == shellNum));
            end
        end
    end
    
end


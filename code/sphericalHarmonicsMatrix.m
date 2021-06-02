classdef sphericalHarmonicsMatrix < matlab.mixin.SetGet
    %SPHERICALHARMONICSMATRIX: this class holds the coefficients and data for conversion to and from SH domain.
    %   Constructor:
    %   sh = sphericalHarmonicsMatrix;
    %   creates an empty object. We need to populate this with the order, compute the T matrix using bvalues and bvecs.
    %   Properties include:
    %   T matrix, pseudo-inverse of the T matrix pT, L, C matrices for the regularization and order to determine the specificlty of the T matrix.
    %   Methods include:
    %   make (to make the T matrix from the bvals, bvecs, and several setters and getters for all these properties above.
    
    properties
        T
        pT
        L
        C
        order
    end
    
    methods
        function SH = sphericalHarmonicsMatrix(gradientDirections, order, lambda)
            %These are the default properties.
            SH.order = order;
            SH.make(gradientDirections, order, lambda);
        end
        
        function set.order(SH, orderIn)
            if(orderIn ~= 4 && orderIn ~= 6 && orderIn ~= 8)
                error('Order has to be either 4, 6 or 8. No other values are allowed!');
            else
                SH.order = orderIn;
            end
        end  
        
        function make(SH, gradientDirections, order, lambda)
            SH.order = order;
            [SH.L, SH.C] = SH.makeL(order);
            SH.createTransformMatrix(gradientDirections, order, lambda);
        end
        
        function data = coeffToSample(SH, coeff)
            % Forward transform through SH.T
            data = SH.T*coeff(:);
        end
        
        function coeff = sampleToCoeff(SH, sample)
            % Inverse transform using pseudoinverse (regularized with
            % lambda).
            coeff = SH.pT*sample(:);
        end
        
        function odf = computeODF(SH, data)
            odf = (diag(SH.C.*SH.L/(16*pi^2))*SH.pT)*data(:);
            odf(1) = (1/4*pi);
        end        
    end
    
    methods (Access = private)
        function [L, C] = makeL(~, order)
            
            L = zeros((order+1)*(order+2)/2, 1); % Zeros of length , column vector.
            C = L;
            
            for k = 0:2:order % 0:2:4, = [0 2 4]
                for m = -k:k % [0], [-2 -1 0 1 2], [-4 -3 -2 -1 0 1 2 3 4], ...
                    j = k*(k+1)/2 + m + 1; % [0], [2 1 2 3/2 2/5 0 0 0 0 0 ...]
                    L(j) = -k*(k+1);
                    if nargout>1
                        C(j) = ((-1)^(k/2))*prod(1:2:(k-1))/prod(2:2:k);
                    end
                end
            end
            C = 2*pi*C;
        end
        
        function createTransformMatrix(SH, gradientDirections, order, lambda)
            %lambda is the regularization parameter that has to be included to correct
            %for the high frequency perturbations in the estimation. If not specified,
            %it is assumed to be 0, and pT is just the pseudoinverse of T.
            
            %T is the same as the matrix B in the Descoteaux paper. This is used for
            %estimation, and reconstruction (using it's pseudoinverse).
            
            nSH = (1/2)*(order+1)*(order+2);
            
            [TH, PHI] = cart2sph(gradientDirections(:,1), gradientDirections(:,2), gradientDirections(:,3));

            angles(:,1) = (pi/2)-PHI;
            angles(:,2) = TH;
            
            nGrad = size(angles,1); % Number of gradient directions.
            SH.T = zeros(nGrad,nSH); % preallocate.
            for j=1:nSH
                v = zeros(1,nSH);
                v(j) = 1;
                for i=1:nGrad
                    % This is the actual SH fitting code. This can be precalculated for a known gradient table.
                    SH.T(i,j) = SH.constructSHBasis(v, angles(i,1), angles(i,2));
                end
            end
            % Calculated as the pseudo inverse of the T matrix.
            SH.pT = inv(SH.T'*SH.T + lambda.*diag(SH.L.^2))*SH.T';
        end
    
        function s = constructSHBasis(~, a, theta, phi)
            L = (sqrt(8*length(a)+1)-3)/2;
            s = 0;
            
            for k=0:2:L
                Pkm=legendre(k, cos(theta))';
                for m=-k:k
                    j=k*(k+1)/2+m+1;
                    if m<0
                        Y=sqrt(((2*k+1)/(2*pi))*factorial(k+m)/factorial(k-m))*Pkm(:,-m+1)*cos(m*phi);
                    elseif m==0
                        Y=sqrt((2*k+1)/(4*pi))*Pkm(:,1)*ones(1,length(phi));
                    else
                        Y=(-1)^m*sqrt(((2*k+1)/(2*pi))*factorial(k-m)/factorial(k+m))*Pkm(:,m+1)*sin(m*phi);
                    end
                    s=s+a(j)*Y;
                end
            end
        end
    end 
end

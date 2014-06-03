classdef SHmat < hgsetget
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T;
        pT;
        L;
        C;
        order;
    end
    
    methods
        function Tmat = SHmat(Tmat)
            Tmat.order = 4;
            Tmat.T = [];
            Tmat.pT = [];
            Tmat.L = [];
            Tmat.C = [];
            
        end
        
        function set.T(Tmat,mat)
            Tmat.T = mat;
        end
        
        function set.pT(Tmat,mat)
            Tmat.pT = mat;
        end
        
        function set.order(Tmat,orderIn)
            if(orderIn ~= 4 && orderIn ~= 6 && orderIn ~= 8)
                error('Order has to be either 4, 6 or 8. No other values are allowed!');
            else
                Tmat.order = orderIn;
            end
        end  
        
        function make(Tmat, bvecs, ord, lambda)
            Tmat.order = ord;
            [Tmat.L,Tmat.C] = Tmat.makeL(ord);             
            Tmat.createPT(bvecs,ord,lambda);           
        end
        
        function data = getData(Tmat, SHCoeff)
            data = Tmat.T*SHCoeff(:);
        end
        
        function coeff = getCoeff(Tmat, data)
             coeff = Tmat.pT*data(:);
        end
        
        function odf = getODF(Tmat, data)
            odf = (diag(Tmat.C.*Tmat.L/(16*pi^2))*Tmat.pT)*data(:);
            odf(1) = (1/4*pi);
        end        
    end
    
    methods (Access = private)
        function [L,C] = makeL(Tmat, ord)
            %Usually order = 4.
            L = zeros((ord+1)*(ord+2)/2, 1); % Zeros of length 15, column vector.
            C = L;
            
            for k = 0:2:ord % 0:2:4, = [0 2 4]
                for m = -k:k % [0], [-2 -1 0 1 2], [-4 -3 -2 -1 0 1 2 3 4]
                    j = k*(k+1)/2 + m + 1; % [0], [2 1 2 3/2 2/5 0 0 0 0 0 ...]
                    L(j) = -k*(k+1);
                    if nargout>1
                        C(j) = ((-1)^(k/2))*prod(1:2:(k-1))/prod(2:2:k);
                    end
                end
            end
            C = 2*pi*C;
        end
        
        function createPT(Tmat, bvecs, ord, lambda)
            %nSH is the number of SH coefficients, equal to 15 if ord = 4.
            %angles is the table for theta and phi, on which the T matrix is
            %calculated.
            %lambda is the regularization parameter that has to be included to correct
            %for the high frequency perturbations in the estimation. If not specified,
            %it is assumed to be 0, and pT is just the pseudoinverse of T.
            
            %T is the same as the matrix B in the Descoteaux paper. This is used for
            %estimation, and reconstruction (using it's pseudoinverse).
            
            %ord = (sqrt(8*nSH+1)-3)/2;
            nSH = (1/2)*(ord+1)*(ord+2);
            
            %angles = (bvecs(:,3)*[1 1 1]).*[sin(bvecs(:,1)).*cos(bvecs(:,2)), sin(bvecs(:,1)).*sin(bvecs(:,2)), cos(bvecs(:,1))];
            [TH, PHI, R] = cart2sph(bvecs(:,1), bvecs(:,2), bvecs(:,3));

            angles(:,1) = pi/2-PHI;
            angles(:,2) = TH;
            angles(:,3) = R;
            
            nGrad = size(angles,1); % Number of gradient directions.
            Tmat.T = zeros(nGrad,nSH); % preallocate.
            for j=1:nSH
                v = zeros(1,nSH); % preallocate v to be 1 x 15 matrix of zeros.
                v(j)=1; % Only that element in v is non-zero. Rest are all zeros.
                for i=1:nGrad
                    % This is the actual SH fitting code. This can be precalculated for a known gradient table.
                    Tmat.T(i,j) = Tmat.SHvol(v, angles(i,1), angles(i,2));
                end
            end
            %pT = pinv(T); % Calculated as the pseudo inverse of the T matrix. 128x15 and 15x128.
            Tmat.pT = inv(Tmat.T'*Tmat.T + lambda.*diag(Tmat.L.^2))*Tmat.T';
        end
    
        function s = SHvol(Tmat, a, theta, phi)
            L=(sqrt(8*length(a)+1)-3)/2;
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
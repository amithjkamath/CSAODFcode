function [dataInt] = IntData(data3X, bvecsI, bvecsO, shellN, order, lambda)

%bvecsI is the bvecs for the input. 
%bvecsO is the bvecs for the output. This can potentially be different.
% when bvecsO is the same as bvecsI, this is effectively interpolating over
% the existing directions, to estimate the samples as an aligned scheme, as
% opposed to a staggered scheme.

nshells = max(shellN);

if ~exist('order','var')
    order = 8; %Default to order 8 reconstruction.
end

if ~exist('lambda','var')
    lambda = zeros(nshells,1);
end

for ns = 1:nshells
    bvecsSH = bvecsI((shellN == ns),:);
    dataSH = data3X(shellN == ns);
    sh = SHmat;
    sh.make(bvecsSH,order,lambda(ns));
    sigSH(ns,:) = sh.getCoeff(dataSH);
end

shO = SHmat;
shO.make(bvecsO,order,0);

for ns = 1:nshells
dataInt(ns,:) = shO.T*sigSH(ns,:)';
end

dataInt = dataInt';

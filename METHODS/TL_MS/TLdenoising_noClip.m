function [Iss_tl,paramsout,runtime] = TLdenoising_noClip(imn,sigma,n,sparsity_controller)
% This function calls "TSPCLOSEDFORMdenoising_noClip"
%
% See "TLdenoising.m" and 
%


% ===================
% Set "C" constant
% ===================
% C: determines sparsity levels in the variable sparsity update step
if nargin<4 || sparsity_controller == -1
    if(n==64)
        C=1.08;
    end
    if(n==121)
        C=1.04;
    end
else
    C = sparsity_controller;
end



sig=sigma;  %standard deviation of Gaussian noise

s=round((0.1)*n);
if(sig>=100)
    iter=5;
else
    iter=11;
end



% initialize parameters

paramsin.sig = sig;
paramsin.iterx = iter;
paramsin.n = n;
paramsin.N = 32000;

paramsin.C = C;
paramsin.tau = 0.01/sig;%0.01/sig
paramsin.s = s;
paramsin.M = 12;%12
paramsin.maxsparsity = round(6*s);
paramsin.method = 0;
paramsin.lambda0 = 0.031;
paramsin.W = kron(dctmtx(sqrt(n)),dctmtx(sqrt(n)));
paramsin.r = 1;



%Denoising algorithm

runtime = cputime;
[Iss_tl,paramsout]= TSPCLOSEDFORMdenoising_noClip(imn,paramsin); %IMU is the output denoised image
runtime = cputime - runtime;


end
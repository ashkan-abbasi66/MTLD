function [Iss_tl,paramsout,runtime] = TLdenoising(imn,im,sigma,n,verbose)
% This function calls the original "TSPCLOSEDFORMdenoising".
%
% INPUTS:
% imn: noisy image
% im: ground-truth image
% sigma: standard deviation of Gaussian noise
% n: patch size. e.g., 64,121,...
% 
% OUTPUTS:
% out.ANh: denoised N-th approximation coefficient
% Iss_tl: denoised image
% paramsout: Structure containing outputs other than the denoised image.
%            paramsout.PSNR : PSNR of denoised image.
%            paramsout.transform : learnt transform
% runtime: elapsed time for the procedure (in seconds)
% 
% USAGE:
%   [Iss_tl,paramsout,runtime] = TLdenoising(imn,im,sigma,64); % fast
%   [Iss_tl,paramsout,runtime] = TLdenoising(imn,im,sigma,121);
% 


if nargin<5
    verbose = 1;
end


if verbose == 1,disp('~~~~ Perform TLD ....'),end

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
if(n==64)
    C=1.08;
end
if(n==121)
    C=1.04;
end
paramsin.C = C;
paramsin.tau = 0.01/sig;
paramsin.s = s;
paramsin.M = 12;
paramsin.maxsparsity = round(6*s);
paramsin.method = 0;
paramsin.lambda0 = 0.031;
paramsin.W = kron(dctmtx(sqrt(n)),dctmtx(sqrt(n)));
paramsin.r = 1;



%Denoising algorithm

runtime = cputime;
[Iss_tl,paramsout]= TSPCLOSEDFORMdenoising(imn,im,paramsin); %IMU is the output denoised image
runtime = cputime - runtime;

end
function [imout, dict] = ...
    ksvd_denoising_Iss(imnoise,sigma,blocksize,iternum,DictSS,lambda)
% This can be used to 
% to obtain "Iss" in the Fused KSVD package
% 

if nargin<5
    iternum = 20;
end

params.x = imnoise;
params.blocksize = blocksize;
params.sigma = sigma;
params.maxval = 255;
params.trainnum = 40000;
params.iternum = iternum;
params.memusage = 'high';

params.initdict = DictSS;
params.lambda = lambda;

msgdelta = -1;


% denoise!
[imout, dict] = ksvddenoise_Iss(params,msgdelta);
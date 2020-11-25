% See "denoise_VST_KSVD13.m"
%
% addpath('ksvdbox13')
% addpath('ompbox10')

% USAGE: (KSVD denoising (ksvdbox13))
%
% blocksize = 8;
% RR=4; % redundancy factor
% dictsize = RR*bb^2; % number of atoms in the dictionary
% IterNum = 20;
% [Iksvd13, ~] = ksvd_denoising(imn,sig,blocksize,dictsize,IterNum);

function [imout, dict] = ksvd_denoising(imnoise,sigma,blocksize,dictsize,iternum)

if nargin<5
    iternum = 20;
end

params.x = imnoise;
params.blocksize = blocksize;
params.dictsize = dictsize;
params.sigma = sigma;
params.maxval = 255;
params.trainnum = 40000;
params.iternum = iternum;
params.memusage = 'high';

msgdelta = -1;


% denoise!
[imout, dict] = ksvddenoise(params,msgdelta);
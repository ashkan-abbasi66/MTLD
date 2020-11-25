% USAGE: (KSVD denoising (2006))
%
% [Iksvd2006, ~] = ksvd_denoising_2006(imn,sigma)

function [Iksvd2006, paramout_ksvd] = ksvd_denoising_2006(imn,sigma)

bb=8; % block size
RR=4; % redundancy factor
K=RR*bb^2; % number of atoms in the dictionary
[Iksvd2006,paramout_ksvd] = denoiseImageKSVD(imn, sigma,K);
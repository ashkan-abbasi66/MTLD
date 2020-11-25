%% Denoising Demo

clear
clc

% Load Image

% load Lena
load ImLandscape   % larger landscape image, where texture artifacts are more noticeable

% Preparing data
NoiseSigma=25;
Imn= Im + NoiseSigma*randn(size(Im));

% loading dictionaries
load DictMS
load DictSS

% Denoising 
[Imd,Iss]=FusedKSVD(Imn,NoiseSigma,DictMS,DictSS,1);

% evaluating
PsnrF = PSNR(Imd,Im); PsnrKsvd = PSNR(Iss,Im);
SsimF = ssim(Imd,Im); SsimKsvd = ssim(Iss,Im);

% Plotting
figure
imshow(Iss,[0 255]);
title(['K-SVD, SSIM = ',num2str(SsimKsvd),', PSNR = ',num2str(PsnrKsvd)]);

figure
imshow(Imd,[0 255]);
title(['Fused K-SVD, SSIM = ',num2str(SsimF),', PSNR = ',num2str(PsnrF)]);

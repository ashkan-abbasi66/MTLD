% INPUTS:
% im: orginal image
% imf: reconstructed image
% imn: noisy image (OPTIONAL)
%
% Usage:
%   PSNR=comp_psnr_2(im,imf)
%   [PSNR,SSIM]=comp_psnr_2(im,imf)
%   [PSNR,SSIM,ISNR]=comp_psnr_2(im,imf,imn)
%
% Ashkan
function [PSNR,SSIM,ISNR]=comp_psnr_2(im,imf,imn,MaxI)

if ~exist('MaxI','var') || isempty(MaxI)
    MaxI = -1;
end
    
if MaxI == -1
    if max(im(:))<2
        MaxI=1;
    else
        MaxI=255;
    end
end

MSE=mean(mean((im(:)-imf(:)).^2));

PSNR=10*log10((MaxI^2)/MSE);
%
if nargout>1
    SSIM=ssim(imf,im,'DynamicRange',MaxI); % [ssimval,ssimmap] = ssim(A,ref)
end
if nargout>2 && (exist('imn','var') || ~isempty(imn))
    ISNR=10*log10((mean(mean((im(:)-imn(:)).^2))/MSE));
end
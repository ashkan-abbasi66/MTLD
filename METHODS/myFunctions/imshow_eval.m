function imshow_eval(im,imf,name,maxI)

    if nargin<4
        maxI = 255;
    end
    
    msg = ', PSNR: %.2fdB, SSIM: %.4f';
    
    
    if isa(imf,'double') == 1 & max(imf(:))<2
        figure,imshow(imf);
    else
        figure,imshow(imf/maxI);
    end
    
    
    msg_ = [name msg];
    [psnr_,ssim_] = comp_psnr(im,imf);
    
    title(sprintf(msg_,psnr_,ssim_));
    disp(sprintf(msg_,psnr_,ssim_))
end
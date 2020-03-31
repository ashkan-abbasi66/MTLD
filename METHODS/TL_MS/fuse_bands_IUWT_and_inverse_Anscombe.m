function [img_denoise] = fuse_bands_IUWT_and_inverse_Anscombe(Iss_tl_out,outs,func)

    if nargin<3 || isempty(func)
        func = @fuse_bands_IUWT; % Astro filterbank
    end

    disp('fuse bands IUWT and inverse Anscombe ...')
    
    scale_shift = outs.scale_shift;
    scale_range = outs.scale_range;
    maxzans = outs.maxzans;
    minzans = outs.minzans;
%     ANh = outs.ANh;
    Ims_tl = outs.Ims_tl;
    a = outs.a;
    sigma = outs.sigma;

    %D = fuse_bands_IUWT_97(Iss_tl_out,Ims_tl,1);
    %D = fuse_bands_IUWT_53(Iss_tl_out,Ims_tl,1);
    %D = fuse_bands_IUWT_79(Iss_tl_out,Ims_tl,1);
    D = func(Iss_tl_out,Ims_tl,1);
    %D = fuse_bands_SWT(Iss_tl_out,Ims_tl,1,'dmey');
    
    D = D/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %%  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse
    
end
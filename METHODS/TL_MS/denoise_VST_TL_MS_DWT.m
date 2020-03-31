function [img_denoise, time,outs] = denoise_VST_TL_MS_DWT(img_raw,sparsity_controller)
% OUTPUTS:
%{
    img_denoise: final output after scaling and applying inverse Anscombe
    transform.

    time

    outs.scale_shift
    outs.scale_range
    outs.maxzans
    outs.minzans
    outs.ANh
    outs.Ims_tl: output before applying scaling and inverse Anscombe tr.
    outs.a
    outs.sigma

%}
    
    if nargin<2
        sparsity_controller = -1;
    end
    
    disp('~~~~ Perform Multi-scale TLD ....')
    
    %% estimate the noise of raw image
    fitparams = estimate_noise(img_raw);
    a = fitparams(1);
    b = fitparams(2);
    if a<0
        a = eps;
    end
    if b<0
        b = eps;
    end
    sigma = sqrt(b);
    
    
    t0=clock;  

    %%  apply forward variance stabilizing transformation
    fz = GenAnscombe_forward(img_raw, sigma, a, 0);
    sigma_den = 1; % Standard-deviation value assumed after variance-stabiliation
    
    % scale the transformed image to [0, 1]
    scale_range = 1;
    scale_shift = (1-scale_range)/2;
    maxzans = max(fz(:));
    minzans = min(fz(:));
    fz = (fz-minzans)/(maxzans-minzans);   sigma_den = sigma_den/(maxzans-minzans);
    fz = fz*scale_range+scale_shift;      sigma_den = sigma_den*scale_range;
    
    %% perform denoising algorithm on transformed images
    
    wname = 'dmey';
    Ns = 1;
    n = 121;
    verbose = 0;
    [out_ms_tl,~] = TLdenoising_MS(fz*255,sigma_den*255,n,Ns,wname,verbose,sparsity_controller);
    Ims_tl = out_ms_tl.Ims_tl;
    ANh = out_ms_tl.ANh;

    
    D = Ims_tl/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %%  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse

    time = etime(clock, t0);
    
    outs.scale_shift = scale_shift;
    outs.scale_range = scale_range;
    outs.maxzans = maxzans;
    outs.minzans = minzans;
    outs.ANh = ANh;
    outs.Ims_tl = Ims_tl;
    outs.a = a;
    outs.sigma = sigma;
    
end
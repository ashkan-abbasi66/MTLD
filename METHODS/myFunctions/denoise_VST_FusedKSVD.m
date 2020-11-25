function [Ifd, Iss , Ims, runtime] = denoise_VST_FusedKSVD(img_raw)
    
    disp('~~~~ KSVD, MS KSVD, and Fused KSVD ....')
    
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
    
    
    %%  apply forward variance stabilizing transformation
    
    start_time = clock;
    
    fz = GenAnscombe_forward(img_raw, sigma, a, 0);
    sigma_den = 1; % Standard-deviation value assumed after variance-stabiliation
    
    % scale the transformed image to [0, 1]
    scale_range = 1;
    scale_shift = (1-scale_range)/2;
    maxzans = max(fz(:));
    minzans = min(fz(:));
    fz = (fz-minzans)/(maxzans-minzans);   sigma_den = sigma_den/(maxzans-minzans);
    fz = fz*scale_range+scale_shift;      sigma_den = sigma_den*scale_range;
    
    time_part1 = etime(clock, start_time);
    
    
    %% perform denoising algorithm on transformed images
    
    if sigma_den>0.15
        weaken = 0.7;% 0.8:31.05,30.94     ;  0.6: 31.35,45  ;0.7: 31.62,58
    else
        weaken = -1;
    end

    
    load DictMS
    load DictSS

%         [Ifd, Iss , Ims]=FusedKSVD(imn,sig,DictMS,DictSS,1);
    verbose = 1;
    [Ifd, Iss , Ims, runtime] = ...
        FusedKSVD_with_runtime(fz*255, sigma_den*255,DictMS,DictSS,verbose,weaken);
    
    % Rescale and compute inverse Anscombe transform
    
    args.scale_shift = scale_shift;
    args.scale_range = scale_range;
    args.maxzans = maxzans;
    args.minzans = minzans;
    args.a = a;
    args.sigma = sigma;
    
    
    start_time = clock;
    
    Iss = rescale_and_inverse(Iss,args);
    
    time_part2 = etime(clock, start_time);
    
    
    Ims = rescale_and_inverse(Ims,args);
    Ifd = rescale_and_inverse(Ifd,args);
    
    runtime.Iss = time_part1 + time_part2 + runtime.Iss;
    % Assuming that "rescale_and_inverse"'s runtime is similar for Iss, 
    %   Ims, and Ifd:
    runtime.Ims = time_part1 + time_part2 + runtime.Ims;
    runtime.Ifd = time_part1 + time_part2 + runtime.Iss + runtime.Ims + runtime.Ifd;
    
end




%% ========================================================================

function img_denoise = rescale_and_inverse(D,args)
    scale_shift = args.scale_shift;
    scale_range = args.scale_range;
    maxzans = args.maxzans;
    minzans = args.minzans;
    a = args.a;
    sigma = args.sigma;
    
    D = D/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse
end
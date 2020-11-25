function [img_denoise, time,outs] = ...
    denoise_VST_MTLD(img_raw,sparsity_controller,Nscale)
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

% Usage patterns:
[img_denoise, time,outs] = denoise_VST_MTLD(img_raw,sparsity_controller,Nscale)
[img_denoise, time,outs] = denoise_VST_MTLD(img_raw,sparsity_controller)
[img_denoise, time,outs] = denoise_VST_MTLD(img_raw,[],Nscale)

%}
    
    if ~exist('Nscale','var') || isempty(Nscale)
        Nscale = 1;
    end
    
    if ~exist('sparsity_controller','var') || isempty(sparsity_controller)
        sparsity_controller = -1;
    end
    
    
    disp('~~~~ MTLD ....')
    
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
        sparsity_controller = sparsity_controller * 0.8;
    end
    
    wname = 'dmey';
    n = 121;
    verbose = 0;
    
    % MTLD
    [out_ms_tl,runtime_ms] = TLD_M(fz*255,sigma_den*255,n,Nscale,wname,verbose,sparsity_controller);
    
    Ims_tl = out_ms_tl.Ims_tl;
    ANh = out_ms_tl.ANh;
    
    start_time = clock;
    
    D = Ims_tl/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %%  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse

    time_part2 = etime(clock, start_time);
    
    
    time.overall = time_part1 + time_part2 + runtime_ms.overall;
    time.ANh = time_part1 + runtime_ms.ANh;
    
    
    outs.scale_shift = scale_shift;
    outs.scale_range = scale_range;
    outs.maxzans = maxzans;
    outs.minzans = minzans;
    outs.ANh = ANh;
    outs.Ims_tl = Ims_tl;
    outs.a = a;
    outs.sigma = sigma;
    
end
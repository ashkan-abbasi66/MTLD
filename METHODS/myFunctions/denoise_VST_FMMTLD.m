function [img_denoise, time] = denoise_VST_FMMTLD(img_raw,sparsity_controller,Nscale,Iss_tl)
% 
% 
% Usage Patterns:
%   [img_denoise, time] = denoise_VST_FMMTLD(img_raw)
%   [img_denoise, time] = denoise_VST_FMMTLD(img_raw,[],Nscale)
%   [img_denoise, time] = denoise_VST_FMMTLD(img_raw,[],[],Iss_tl)
%   [img_denoise, time] = denoise_VST_FMMTLD(img_raw,[],[],-1)
%   [img_denoise, time] = denoise_VST_FMMTLD(img_raw,sparsity_controller,Nscale,Iss_tl)
% 
% 
    if ~exist('Nscale','var') || isempty(Nscale)
        Nscale = 1;
    end
    
    if ~exist('sparsity_controller','var') || isempty(sparsity_controller)
        sparsity_controller = -1;
    end
    
    if ~exist('Iss_tl','var') || isempty(Iss_tl)
        Iss_tl = -1;
    end
    
    %disp('~~~~ FMMTLD ....')
    
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
    
    
    % FMMTLD 
    
    n = 121;
    [tl_direct,runtime_fmm] = ...
        TLD_FMM(fz*255,sigma_den*255,n,...
        Nscale,Iss_tl,'',sparsity_controller);  
    
    start_time = clock;
    
    D = tl_direct/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %%  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse
    
    time_part2 = etime(clock, start_time);
    
    time = time_part1 + time_part2 + runtime_fmm;
    
end
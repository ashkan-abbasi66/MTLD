function [Ims, runtime] = denoise_VST_BLS_GSM(img_raw)
    
    disp('~~~ BLS-GSM ....')
    
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
        weaken = 0.7;%0.8;
    else
        weaken = -1;
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    [Ny,Nx] = size(fz);

    PS = ones(size(fz));	% power spectral density (in this case, flat, i.e., white noise)
    seed = 0;               % random seed

    % Pyramidal representation parameters
    Nsc = ceil(log2(min(Ny,Nx)) - 4);  % Number of scales (adapted to the image size)
    Nor = 3;				            % Number of orientations (for X-Y separable wavelets it can only be 3)
    repres1 = 'uw';                     % Type of pyramid (shift-invariant version of an orthogonal wavelet, in this case)
    repres2 = 'daub1';                  % Type of wavelet (daubechies wavelet, order 2N, for 'daubN'; in this case, 'Haar')

    % Model parameters (optimized: do not change them unless you are an advanced user with a deep understanding of the theory)
    blSize = [3 3];	    % n x n coefficient neighborhood of spatial neighbors within the same subband
                        % (n must be odd): 
    parent = 0;			% including or not (1/0) in the neighborhood a coefficient from the same spatial location
                        % and orientation as the central coefficient of the n x n neighborhood, but
                        % next coarser scale. Many times helps, but not always.
    boundary = 1;		% Boundary mirror extension, to avoid boundary artifacts 
    covariance = 1;     % Full covariance matrix (1) or only diagonal elements (0).
    optim = 1;          % Bayes Least Squares solution (1), or MAP-Wiener solution in two steps (0)

    % -------------------------
    % Uncomment the following 4 code lines for reproducing the results of our IEEE Trans. on Im. Proc., Nov. 2003 paper
    % This configuration is slower than the previous one, but it gives slightly better results (SNR)
    % on average for the test images "lena", "barbara", and "boats" used in the cited article.

    Nor = 8;                           % 8 orientations
    repres1 = 'fs';                    % Full Steerable Pyramid, 5 scales for 512x512
    repres2 = '';                      % Dummy parameter when using repres1 = 'fs'   
    parent = 1;                        % Include a parent in the neighborhood
    % -------------------------


    % Call the denoising function
    tic
    if weaken == -1
        Ims = denoi_BLS_GSM(fz*255, sigma_den*255, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);
    else
        Ims = denoi_BLS_GSM(fz*255, sigma_den*weaken*255, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);
    end
    den_time = toc;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Rescale and compute inverse Anscombe transform
    
    args.scale_shift = scale_shift;
    args.scale_range = scale_range;
    args.maxzans = maxzans;
    args.minzans = minzans;
    args.a = a;
    args.sigma = sigma;
    
    
    start_time = clock;
    
    
    time_part2 = etime(clock, start_time);
    
    
    Ims = rescale_and_inverse(Ims,args);
    
    % Assuming that "rescale_and_inverse"'s runtime is similar for Iss, 
    %   Ims, and Ifd:
    runtime = time_part1 + time_part2 + den_time;
    
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
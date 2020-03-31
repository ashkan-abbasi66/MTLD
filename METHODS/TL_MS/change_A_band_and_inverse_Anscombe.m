function [img_denoise] = change_A_band_and_inverse_Anscombe(Iss_tl,outs)

    disp('change A band and inverse Anscombe ...')
    
    scale_shift = outs.scale_shift;
    scale_range = outs.scale_range;
    maxzans = outs.maxzans;
    minzans = outs.minzans;
    ANh = outs.ANh;
    % Ims_tl = outs.Ims_tl;
    a = outs.a;
    sigma = outs.sigma;

    
    wname = 'dmey';
    Ns = 1;
    [D,~] = change_A_band_DWT(Iss_tl,ANh,Ns,wname,0); 
    
    D = D/255;
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;
    
    %%  apply the inverse VST transform
    img_denoise = GenAnscombe_inverse_exact_unbiased(D, sigma, a, 0);   % exact unbiased inverse
    
end
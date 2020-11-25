function Imx = fuse_bands_IUWT_old(Iss,Ims,Ns,func)
% Fuse bands with IUWT (Isotropic Undecimated Wavelet Transform)
% 
% INPUTS:
% Ns: number of decomposition scales. E.g., 1,2, ... 
% wname:  the wavelet name

    if nargin<4
        func = @iuwt;
    end
    
    % Decomposition and mixing

    levels=1:Ns;

    [w_Iss,~]=func(Iss,levels);
    [~,s_out_Ims]=func(Ims,levels);


    % Recomposition

    details = 0;
    for i = 1:Ns
        details = details + w_Iss(:,:,i);
    end

    Imx = s_out_Ims + details;
    
    
end
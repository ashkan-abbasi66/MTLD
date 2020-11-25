function Imx = fuse_bands_IUWT(Iss,Ims,Ns,func)
% Fuse bands with IUWT (Isotropic Undecimated Wavelet Transform)
% 
% INPUTS:
% Ns: number of decomposition scales. E.g., 1,2, ... 
% wname:  the wavelet name

    
    if ~exist('Ns','var') || isempty(Ns)
        Ns = getNs(Iss);
    end
    
    if ~exist('func','var') || isempty(func)
        func = @iuwt;
    end
    
    
    Imx = Iss;
    
    for scale_number = Ns:-1:1

        % Decomposition and mixing

        levels=1:scale_number;

        [w_Iss,~]=func(Imx,levels);
        [~,s_out_Ims]=func(Ims,levels);


        % Recomposition

        details = 0;
        for i = 1:scale_number
            details = details + w_Iss(:,:,i);
        end

        Imx = s_out_Ims + details;
    end
    

end
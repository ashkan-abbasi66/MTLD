function Imx = fuse_bands_IUWT_79(Iss,Ims,Ns)
% Fuse bands with IUWT (Isotropic Undecimated Wavelet Transform)
% 
% INPUTS:
% Ns: number of decomposition scales. E.g., 1,2, ... 
% wname:  the wavelet name



% Decomposition and mixing

levels=1:Ns;

[w_Iss,s_out_Iss]=iuwt_79(Iss,levels);
[w_Ims,s_out_Ims]=iuwt_79(Ims,levels);


% Recomposition

details = 0;
for i = 1:Ns
    details = details + w_Iss(:,:,i);
end

Imx = s_out_Ims + details;


end
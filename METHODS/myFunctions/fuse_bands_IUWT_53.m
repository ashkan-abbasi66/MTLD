function Imx = fuse_bands_IUWT_53(Iss,Ims,Ns)
% Fuse bands with IUWT (Isotropic Undecimated Wavelet Transform)
% 
% INPUTS:
% Ns: number of decomposition scales. E.g., 1,2, ... 
% wname:  the wavelet name


    Imx = fuse_bands_IUWT(Iss,Ims,Ns,@iuwt_53);

end
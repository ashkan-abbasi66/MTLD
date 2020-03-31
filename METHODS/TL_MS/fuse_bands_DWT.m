function Imx = fuse_bands_DWT(Iss,Ims,Ns,wname)
% Fuse bands with DWT (discrete 2-D wavelet transform)
% 
% INPUTS:
% Ns: number of decomposition scales. E.g., 1,2, ... 
%   ("wmaxlev" gives the maximum Ns)
% wname:  the wavelet name



% Decomposition and mixing

[C_ss,S_ss] = wavedec2(Iss,Ns,wname); 
[C_ms,~] = wavedec2(Ims,Ns,wname); 

LL_indices = 1:prod(S_ss(1,:));

C_ss(1,LL_indices) = C_ms(1,LL_indices);

% Recomposition

Imx = waverec2(C_ss,S_ss,wname);


end

% % % % USING DWT2 FUNCTION FOR A 2 LEVEL DECOMPOSITION %
% % % [Iss_cA_1,Iss_cH_1,Iss_cV_1,Iss_cD_1] = dwt2(Iss,wname);
% % % [~,Iss_cH_2,Iss_cV_2,Iss_cD_2] = dwt2(Iss_cA_1,wname);
% % % 
% % % [Ims_cA_1,~,~,~] = dwt2(Ims,wname);
% % % [Ims_cA_2,~,~,~] = dwt2(Ims_cA_1,wname);
% % % 
% % % sz = size(Iss);
% % % imout_level_1 = idwt2(Ims_cA_2,Iss_cH_2,Iss_cV_2,Iss_cD_2,wname,sz);
% % % Imx = idwt2(imout_level_1,Iss_cH_1,Iss_cV_1,Iss_cD_1,wname,sz);
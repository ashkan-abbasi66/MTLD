function [Imx_dwt,runtime] = ...
    change_A_band_DWT(Iss_tl,A1h,Ns,wname,verbose)
% Change approximation coefficient (A1) in DWT
%
% Given one input image (Iss_tl) and one approximation coefficient(A1h), 
%  this function changes the 1st approximation coefficient of Iss_tl
%  with the given approximation coefficient. The decomposition is based
%  on Discrete Wavelet Transform.
%
% INPUTS:
%   Ns: number of decomposition scales
%
% OUTPUTS:
% Imx_tl: the output mixed band image
% runtime: total time of the whole algorithm.
%
% USAGE:
%   [Imx_dwt,runtime] = change_A_band_DWT(Iss_tl_org,out.A1h,1,'db5',1);


if nargin < 5 
    verbose = 1;
end

if verbose == 1,disp('   Change Approximation band in DWT ...'),end

runtime = cputime;

[C_ss,S_ss] = wavedec2(Iss_tl,Ns,wname); 
LL_indices = 1:prod(S_ss(1,:));
C_ss(1,LL_indices) = A1h(:)';

Imx_dwt = waverec2(C_ss,S_ss,wname);

runtime = cputime - runtime;
function [out,runtime] = TLD_M(imn,sigma,n,varargin)
% 
% MTLD (Multi-scale TL denoising)
%
% INPUTS:
% imn: noisy image
% sigma: standard deviation of Gaussian noise
% n: patch size. e.g., 64,121,...
% Ns: number of decomposition scales. E.g., 1,2, ...
% wname: wavelet name. e.g., 'db5'
% verbose: set to 1 for displaying messages
%
% OUTPUTS:
% out.ANh: denoised N-th approximation coefficient
% out.Ims_tl: denoised image using multi-scale transform learning.
% runtime.ANh: elapsed time for DWT and denoising N-th approx. coefficient.
% runtime.overall: elapsed time for the procedure (in seconds).
%
% USAGE:
% n = 121;
% Ns = 2;
% wname = 'db5';
% [out,runtime] = TLD_M(imn,sigma,n,Ns,wname,1);
% 
% Some Notes and Usage Patterns:
%   (1) [out,runtime] = TLD_M(imn,sigma,n)
%   (2) [out,runtime] = TLD_M(imn,sigma,n,Ns)
%   (3) [out,runtime] = TLD_M(imn,sigma,n,Ns,wname,verbose,sparsity_controller);
%   wname can be set to ''
%   verbose can be set to -1 or 1
%   If you do not set sparsity_controller, it is set to its default value
% 
% Ashkan 


% Parse input arguments

p = inputParser;
% validScalarN = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(p,'imn');
addRequired(p,'sigma');
addRequired(p,'n');

addOptional(p,'Ns',getNs(imn)) % (...,validScalarN)
addOptional(p,'wname','dmey',@ischar)
addOptional(p,'verbose',1,@isscalar);
addOptional(p,'sparsity_controller',-1,@isscalar);


parse(p,imn,sigma,n,varargin{:});

Nscale = p.Results.Ns;
wname = p.Results.wname;
verbose = p.Results.verbose;
sparsity_controller = p.Results.sparsity_controller;

if isempty(wname)
    wname = 'dmey';
end

if isempty(Nscale) || Nscale <0
    Nscale = getNs(imn);
end


msg = '~~~~ MTLD ....';
if verbose == 1,disp(msg),end



% Initialize outputs
out.ANh = 0;
out.Ims_tl = 0;
runtime.ANh = 0;
runtime.overall = 0;



% Decomposition

if verbose == 1
    fprintf('   perform DWT with %s and Ns = %d level(s) ...\n',wname,Nscale)
end

tic
[C,L]=wavedec2(imn,Nscale,wname);
dwt_time = toc;

Crec=[];  
detail_times = 0;




% Denoising coefficients

for s=Nscale:-1:1
    
    sig_ = sigma;
    
    if s==Nscale
        msg = '   run TLD for "A%d" ...\n';
        if verbose == 1,fprintf(msg,s),end
        
        AN = appcoef2(C,L,wname,s);
        tic
        [out.ANh,~] = TLD_2_noClip(AN,sig_,n,sparsity_controller);
        runtime.ANh = toc;
        
        Crec=[Crec,out.ANh(:)'];        
    end
    
    
    [H,V,D] = detcoef2('all',C,L,s);
    
    
    msg = '   run TLD for "H%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    
    tic
    [H_hat,~] = TLD_2_noClip(H,sig_,n,sparsity_controller);
    H_time = toc;
    
    
    msg = '   run TLD for "V%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    
    tic
    [V_hat,~] = TLD_2_noClip(V,sig_,n,sparsity_controller);
    V_time = toc;
    
    
    msg = '   run TLD for "D%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    
    tic
    [D_hat,~] = TLD_2_noClip(D,sig_,n,sparsity_controller);
    D_time = toc;
    

    Crec=[Crec,H_hat(:)',V_hat(:)',D_hat(:)'];
    
    detail_times = detail_times + H_time + V_time + D_time;
    
end



if verbose == 1,disp('   inverse DWT ...'),end

tic
out.Ims_tl=waverec2(Crec,L,wname);
runtime.overall = toc + detail_times + runtime.ANh;

runtime.ANh = runtime.ANh + dwt_time;

end



% =================================================
% =================================================

function Ns = getNs(imn)
    % Get the number of decomposition levels
    
    [R,C] = size(imn);
    lmin = min(R,C);
    if lmin>550
        if lmin > 1200
            Ns = 3;
        else
            Ns=2;
        end
    else
        Ns=1;
    end

end
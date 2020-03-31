function [out,runtime] = TLdenoising_MS(imn,sigma,n,varargin)
% Multi-scale TL denoising
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
% [out,runtime] = TLdenoising_MS(imn,sigma,n,Ns,wname,1);
% 
% Ashkan 


% Parse input arguments

p = inputParser;
validScalarN = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(p,'imn');
addRequired(p,'sigma');
addRequired(p,'n');
addOptional(p,'Ns',getNs(imn),validScalarN)
addOptional(p,'wname','dmey',@ischar)
addOptional(p,'verbose',1,@isscalar);
addOptional(p,'sparsity_controller',-1,@isscalar);


parse(p,imn,sigma,n,varargin{:});

verbose = p.Results.verbose;
wname = p.Results.wname;
Ns = p.Results.Ns;
sparsity_controller = p.Results.sparsity_controller;


msg = '~~~~ Perform Multi-scale TLD ....';
if verbose == 1,disp(msg),end



% Initialize outputs
out.ANh = 0;
out.Ims_tl = 0;
runtime.ANh = 0;
runtime.overall = 0;



% Decomposition

if verbose == 1
    fprintf('   perform DWT with %s and Ns = %d level(s) ...\n',wname,Ns)
end

dwt_time = cputime;
[C,L]=wavedec2(imn,Ns,wname);
dwt_time = cputime - dwt_time;

Crec=[];  
overall_time = 0;




% Denoising coefficients

for s=Ns:-1:1
    
    
    if s==Ns
        msg = '   run TLD for "A%d" ...\n';
        if verbose == 1,fprintf(msg,s),end
        
        AN = appcoef2(C,L,wname,s);
        A_time = cputime;
        [out.ANh,~] = TLdenoising_noClip(AN,sigma,n,sparsity_controller);
        runtime.ANh = cputime - A_time;
        runtime.ANh = runtime.ANh + dwt_time;
        
        Crec=[Crec,out.ANh(:)'];        
    end
    
    
    [H,V,D] = detcoef2('all',C,L,s);
    
    
    msg = '   run TLD for "H%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    
    H_time = cputime;
    [H_hat,~] = TLdenoising_noClip(H,sigma,n,sparsity_controller);
    H_time = cputime - H_time;
    
    
    msg = '   run TLD for "V%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    V_time = cputime;
    [V_hat,~] = TLdenoising_noClip(V,sigma,n,sparsity_controller);
    V_time = cputime - V_time;
    
    
    msg = '   run TLD for "D%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    D_time = cputime;
    [D_hat,~] = TLdenoising_noClip(D,sigma,n,sparsity_controller);
    D_time = cputime - D_time;
    

    Crec=[Crec,H_hat(:)',V_hat(:)',D_hat(:)'];
    
    overall_time = overall_time + H_time + V_time + D_time;
    
end

overall_time = overall_time + runtime.ANh;


if verbose == 1,disp('   inverse DWT ...'),end

tic
out.Ims_tl=waverec2(Crec,L,wname);
runtime.overall = overall_time + toc;

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
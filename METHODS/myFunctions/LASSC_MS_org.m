% 
% addpath('E:\MATLAB_CODES\MATLAB\TL_MTLD_SeeFdrive\METHODS\LASSC_Denoising_2')
% 
function [out,runtime] = LASSC_MS_org(imn,sigma,varargin)
% Multi-scale LASSC denoising 
%
% INPUTS:
% imn: noisy image
% sigma: standard deviation of Gaussian noise

% wname: wavelet name. e.g., 'db5'
% verbose: set to 1 for displaying messages
%
% OUTPUTS:
% out.ANh: denoised N-th approximation coefficient
% out.im_out : denoised image using multi-scale transform learning.
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
addOptional(p,'Ns',getNs(imn),validScalarN)
addOptional(p,'wname','dmey',@ischar)
addOptional(p,'verbose',1,@isscalar);
% addOptional(p,'sparsity_controller',-1,@isscalar);


parse(p,imn,sigma,varargin{:});

verbose = p.Results.verbose;
wname = p.Results.wname;
Ns = p.Results.Ns;
% sparsity_controller = p.Results.sparsity_controller;


msg = '============ Multi-scale LASSC ============';
if verbose == 1,disp(msg),end



% Initialize outputs
out.ANh = 0;
out.im_out  = 0;
runtime.ANh = 0;
runtime.overall = 0;



% Decomposition

if verbose == 1
    fprintf('   DWT with %s and Ns = %d level(s) ...\n',wname,Ns)
end

tic
[C,L]=wavedec2(imn,Ns,wname);
dwt_time = toc;

Crec=[];  
overall_time = 0;




% Denoising coefficients

for s=Ns:-1:1
    
    
    if s==Ns
        msg = '   LASSC denoising for "A%d" ...\n';
        if verbose == 1,fprintf(msg,s),end
        
        AN = appcoef2(C,L,wname,s);
        tic
%         [out.ANh,~] = TLdenoising_noClip(AN,sigma,n,sparsity_controller);
        out.ANh = run_LASSC_denoising_org(AN,sigma);
        runtime.ANh = toc;
        runtime.ANh = runtime.ANh + dwt_time;
        
        Crec=[Crec,out.ANh(:)'];        
    end
    
    
    [H,V,D] = detcoef2('all',C,L,s);
    
    
    msg = '   LASSC denoising for "H%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    
    tic
%     [H_hat,~] = TLdenoising_noClip(H,sigma,n,sparsity_controller);
    H_hat = run_LASSC_denoising_org(H,sigma);
    H_time = toc;
    
    
    msg = '   LASSC denoising for "V%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    tic
%     [V_hat,~] = TLdenoising_noClip(V,sigma,n,sparsity_controller);
    V_hat = run_LASSC_denoising_org(V,sigma);
    V_time = toc;
    
    
    msg = '   LASSC denoising for "D%d" ...\n';
    if verbose == 1,fprintf(msg,s),end
    tic
%     [D_hat,~] = TLdenoising_noClip(D,sigma,n,sparsity_controller);
    D_hat = run_LASSC_denoising_org(D,sigma);
    D_time = toc;
    

    Crec=[Crec,H_hat(:)',V_hat(:)',D_hat(:)'];
    
    overall_time = overall_time + H_time + V_time + D_time;
    
end

overall_time = overall_time + runtime.ANh;


if verbose == 1,disp('   Inverse DWT ...'),end

tic
out.im_out =waverec2(Crec,L,wname);

runtime.overall = toc + overall_time;

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
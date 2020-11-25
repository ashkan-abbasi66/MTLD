function [out,runtime] = ...
    TLD_FMM(imn,sigma,n,varargin)
%   
%   FMMTLD
% 
% runtime: note that it does not consider the runtime of single-scale image
% denoising when "Iss_tl" is given as input.
% 
% USAGE:
%   (1)
%   Nscale = 2;
%   out = TLD_FMM(imn,sig,121,Nscale);
% 
%   If you want to just give one of the optional arguments, 
%      use these commands:
% 
%   (2)
%   out = TLD_FMM(imn,sig,121,Nscale,Iss_tl);
% 
%   (3)
%   out = TLD_FMM(imn,sig,121,Nscale,-1,'dmey');
%   
%   (4)
%   sparsity_controller = 1.04;
%   out = TLD_FMM(imn,sig,121,Nscale,-1,'',sparsity_controller);
% 
% 

total_time = 0;

disp('~~~~ FMMTLD ...');

% Parse input arguments

p = inputParser;
validScalarN = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(p,'imn');
addRequired(p,'sigma');
addRequired(p,'n');

addOptional(p,'Nscale',getNs(imn),validScalarN)
addOptional(p,'Iss_tl',-1,@isnumeric);
addOptional(p,'wname','dmey',@ischar)
addOptional(p,'sparsity_controller',-1,@isscalar);


parse(p,imn,sigma,n,varargin{:}); %%%%

Nscale = p.Results.Nscale;
Iss_tl = p.Results.Iss_tl;
wname = p.Results.wname;
sparsity_controller = p.Results.sparsity_controller;

if isempty(wname)
    wname = 'dmey';
end


% =================================

tic

LL{1} = imn;
for i = 1:Nscale
    [A,~,~,~] = dwt2(LL{i},wname);
    LL{i+1} = A;
end

total_time = total_time + toc;


% % % DEBUG
% % % for i = 1:Nscale+1
% % %     disp(i)
% % %     sigma_hat = function_stdEst(LL{i})
% % % end



% Denoising coefficients
LL_hat = LL;


% Denoise Last LL coefficient
ind = Nscale + 1;
fprintf('    Denoising LL%d\n',ind)
[LL_hat{ind},~,rt] = TLD_2_noClip(LL{ind},sigma,n,sparsity_controller);

total_time = total_time + rt;

% Denoise remaining LL coefficient
for i = Nscale:-1:1
    
    ind = i+1;
    tmpA = LL_hat{ind};
    
    ind = i;
    fprintf('    Denoising LL%d\n',ind)
    if ind == 1        
        if Iss_tl ~= -1 
            % When the result of single-scale denoising is given by user
            LL_hat{ind} = Iss_tl;
            rt = 0; % You should add it later.
        else
            [LL_hat{ind},~,rt] = TLD_2(LL{ind},sigma,n,sparsity_controller); 
        end
        
    else
        [LL_hat{ind},~,rt] = TLD_2_noClip(LL{ind},sigma,n,sparsity_controller); 
    end
    
    
    
    tic
    
    % fusion
    [~,H,V,D] = dwt2(LL_hat{ind},wname);
    tmpA = tmpA(1:size(H,1),1:size(H,2));
    LL_hat{ind} = idwt2(tmpA,H,V,D,wname);
    
    total_time = toc + total_time + rt;
end



out = LL_hat{1};

runtime = total_time;

end




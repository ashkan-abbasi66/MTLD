function [out] = ...
    WNNM_FMM(imn,sigma,varargin)

disp('~~~~ WNNM - FMM ...');

% Parse input arguments

p = inputParser;
validScalarN = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(p,'imn');
addRequired(p,'sigma');

addOptional(p,'Nscale',getNs(imn),validScalarN)
addOptional(p,'Iss_tl',-1,@isnumeric);
addOptional(p,'wname','dmey',@ischar)
addOptional(p,'sparsity_controller',-1,@isscalar);


parse(p,imn,sigma,varargin{:}); %%%%

Nscale = p.Results.Nscale;
Iss_tl = p.Results.Iss_tl;
wname = p.Results.wname;
sparsity_controller = p.Results.sparsity_controller;

if isempty(wname)
    wname = 'dmey';
end


% =================================


LL{1} = imn;
for i = 1:Nscale
    [A,~,~,~] = dwt2(LL{i},wname);
    LL{i+1} = A;
end

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
% [LL_hat{ind},~,runtime] = TLD_2_noClip(LL{ind},sigma,n,sparsity_controller);
Par = ParSet(sigma);
LL_hat{ind} = WNNM_DeNoising(LL{ind}, zeros(size(LL{ind})), Par); 



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
            
        else
            %[LL_hat{ind},~,runtime] = TLD_2(LL{ind},sigma,n,sparsity_controller); 
			%Par = ParSet(sigma);
			LL_hat{ind} = WNNM_DeNoising(LL{ind}, zeros(size(LL{ind})), Par); 
        end
        
    else
        %[LL_hat{ind},~,runtime] = TLD_2_noClip(LL{ind},sigma,n,sparsity_controller); 
		LL_hat{ind} = WNNM_DeNoising(LL{ind}, zeros(size(LL{ind})), Par);
    end
    
    
    
    
    
    % fusion
    [~,H,V,D] = dwt2(LL_hat{ind},wname);
    tmpA = tmpA(1:size(H,1),1:size(H,2));
    LL_hat{ind} = idwt2(tmpA,H,V,D,wname);
    
    
end



out = LL_hat{1};



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



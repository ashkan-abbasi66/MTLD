function out = TLdenoising_MS_dwtFused(imn,wname,Nscale,sigma)%,sigma,n,wname,verbose)


% addpath('E:\mfiles_acode_thesis\SBSDI_2013_Fang_2');

LL{1} = imn;
for i = 1:Nscale
    [A,~,~,~] = dwt2(LL{i},wname);
    LL{i+1} = A;
end


for i = 1:Nscale+1
    disp(i)
    sigma_hat = function_stdEst(LL{i})
end



% Denoising coefficients
LL_hat = LL;

n = 121;


% Denoise Last LL coefficient
ind = Nscale + 1;
fprintf('    Denoising LL%d\n',ind)
[LL_hat{ind},~] = TLdenoising_noClip(LL{ind},sigma,n);


% Denoise remaining LL coefficient
for i = Nscale:-1:1
    
    ind = i+1;
    tmpA = LL_hat{ind};
    
    ind = i;
    fprintf('    Denoising LL%d\n',ind)
    [LL_hat{ind},~] = TLdenoising_noClip(LL{ind},sigma,n); 
    
    % fusion
    [~,H,V,D] = dwt2(LL_hat{ind},wname);
    LL_hat{ind} = idwt2(tmpA,H,V,D,wname);
end



out = LL_hat{1};







% % % 
% % % if verbose == 1,disp('   TL Denoising of "A1" ...'),end
% % % 
% % % A1h_time = cputime;
% % % [out.A1h,~] = TLdenoising_noClip(A1,sigma,n);
% % % runtime.A1h = cputime - A1h_time;
% % % runtime.A1h = runtime.A1h + dwt_time;
% % % 
% % % 
% % % if verbose == 1,disp('   TL Denoising of "H1" ...'),end
% % % 
% % % H1h_time = cputime;
% % % [H1h,~] = TLdenoising_noClip(H1,sigma,n);
% % % H1h_time = cputime - H1h_time;
% % % 
% % % 
% % % if verbose == 1,disp('   TL Denoising of "V1" ...'),end
% % % 
% % % V1h_time = cputime;
% % % [V1h,~] = TLdenoising_noClip(V1,sigma,n);
% % % V1h_time = cputime - V1h_time;
% % % 
% % % 
% % % if verbose == 1,disp('   TL Denoising of "D1" ...'),end
% % % 
% % % D1h_time = cputime;
% % % [D1h,~] = TLdenoising_noClip(D1,sigma,n);
% % % D1h_time = cputime - D1h_time;
% % % 
% % % 
% % % 
% % % % Recomposition
% % % 
% % % if verbose == 1,disp('   Inverse DWT ...'),end
% % % 
% % % idwt_time = cputime;
% % % sz = size(imn);
% % % out.Ims_tl = idwt2(out.A1h,H1h,V1h,D1h,wname,sz);
% % % idwt_time = cputime - idwt_time;
% % % 
% % % 
% % % runtime.overall = runtime.A1h + ... % includes time for performing DWT
% % %     H1h_time + V1h_time + D1h_time + idwt_time;

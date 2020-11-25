function [Iss , Ims, Imx_dwt, Imx_iuwt, runtime ] = MMKSVD_with_runtime(Imn,NoiseSigma,DictMS,DictSS,verbose,weaken,Nscale)
%
% runtime.Imx_dwt
% runtime.Imx_iuwt
% 
%
% Usage Example:
%   See "Benchmark_MMKSVD_for_Gaussian_denoising"
% 
    
    if ( exist('omp2','file')==0 ||  exist('ksvd','file')==0 )
        error('ksvd and omp packages needed.')
    end
    
    if ~exist('weaken','var') || isempty(weaken)
        weaken = -1;
    end
    
    if ~exist('Nscale','var') || isempty(Nscale)
        S = getNscale(Imn);
    else
        S = Nscale;
    end


    
    % Wavelet name:
    Wname = 'dmey';
    
    % K-SVD iterations
    IterNum=20;
   
    [Iss, Ims , runtime]=MScaleKSVDdenoising(S,IterNum,DictMS,DictSS,Wname,Imn,NoiseSigma,verbose,weaken);              
    
    if verbose == 1
        disp('    Fusing...')
    end
    
    
%     n = 121;
%     [Imx_dwt,runtime_fmm] = ...
%         TLD_FMM(imn,sig,n,mx_Ns,Iss_tl,mx_wname); 
    tic
    Imx_dwt = fuse_bands_DWT(Iss,Ims,S,Wname);
    runtime.Imx_dwt = toc + runtime.Iss + runtime.approx_denoising_time;


    tic
    Imx_iuwt = fuse_bands_IUWT(Iss,Ims,S);
    runtime.Imx_iuwt = toc + runtime.Iss + runtime.Ims;
    
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Iss, Ims , runtime] = ...
    MScaleKSVDdenoising(S,IterNum,DictMS,DictSS,Wname,image,NoiseSigma,verbose,weaken_factor)
% Multi scale and single scale K-SVD dict adaptive image denoising by K-SVD
% 
% 
% runtime.Iss = toc;
% runtime.Ims
% runtime.approx_denoising_time
% 

TS = dwtmode('status','nodisp');
if ~strcmp(TS,'sym')
    dwtmode('sym','nodisp')
end

[n,m] = size(DictSS);

if weaken_factor == -1
    E = sqrt(NoiseSigma^2*n);   % coding threshold
else
    E = sqrt(NoiseSigma^2*n)*weaken_factor;   % coding threshold
end

tic

[C,L]=wavedec2(image,S,Wname);

Crec=[];  

if verbose==1
    disp('    Denoising Regular K-SVD...')
end

[Iss]=CleanKSVD(image,n,m,E*getConstant(n),IterNum,DictSS,20/NoiseSigma);
runtime.Iss = toc;

tic

% Preparing MS Dict %

scales = length(DictMS);

DictMs.H = [];
DictMs.D  = [];
DictMs.V  = [];
DictMs.A  = [];

for s=1:scales
    DictMs.H = [DictMs.H, DictMS(s).H];
    DictMs.V = [DictMs.V, DictMS(s).V];
    DictMs.D = [DictMs.D, DictMS(s).D];    
end

DictMs.A = DictMS(scales).A;

% MS Denoising %
if verbose==1
    disp('    Denoising Multi-scale K-SVD...')
end

approx_denoising_time = 0;

for s=S:-1:1
    
    if s==S
        tic
        A = appcoef2(C,L,Wname,s);
        [Ac,~,~] = CleanKSVD(A,n,m,E*getConstant(n*4*s^2,0.95),IterNum,DictMs.A,20/NoiseSigma);
        Crec=[Crec,Ac(:)'];        
        approx_denoising_time = toc + approx_denoising_time;
    end
    
    [H,V,D] = detcoef2('all',C,L,s);
    
    [Vc,~,~] = CleanKSVD(V,n,m,E*getConstant(n*4*s^20,0.999),IterNum,DictMs.V,0);
    
    [Hc,~,~] = CleanKSVD(H,n,m,E*getConstant(n*4*s^2,0.999),IterNum,DictMs.H,0);
    
    [Dc,~,~] = CleanKSVD(D,n,m,E*getConstant(n*4*s^2,0.999),IterNum,DictMs.D,0);
    
    Crec=[Crec,Hc(:)',Vc(:)',Dc(:)'];
        
end

Ims=waverec2(Crec,L,Wname);

runtime.Ims = toc;
runtime.approx_denoising_time = approx_denoising_time;

end


% ========================================================
% ========================================================
% ========================================================

function [Aout,Dict,X,Yd] = CleanKSVD(A,n,m,E,IterNum,InitDict,lamda)

[N,M]=size(A);
p=sqrt(n);
Y=zeros(n,(N-p+1)*(M-p+1));
pt=1;
for i=1:1:N-p+1
    for j=1:1:M-p+1
        patch=A(i:i+p-1,j:j+p-1);
        Y(:,pt)=patch(:);
        pt=pt+1;
    end
end

if isempty(InitDict)
    params.initdict=odct2dict([sqrt(n),sqrt(n)],[round(sqrt(m)),round(sqrt(m))]);
else
    params.initdict=InitDict;
end
params.Edata=E;
params.iternum=IterNum;
params.data=Y;

[Dict,X]=ksvd(params,'');

Yd=Dict*X;

i=1;
j=1;
Weight=zeros(size(A));
Aout=zeros(size(A));
for k=1:size(Yd,2)
    patch=reshape(Yd(:,k),[p,p]);
    Aout(i:i+p-1,j:j+p-1)=Aout(i:i+p-1,j:j+p-1)+patch;
    Weight(i:i+p-1,j:j+p-1)=Weight(i:i+p-1,j:j+p-1)+1;
    if j<M-p+1 
        j=j+1; 
    else
        j=1; i=i+1; 
    end;
end;    
Aout=(lamda*A+Aout)./(lamda+Weight);

end


% ==============================================================
% ==============================================================
% ==============================================================
function S = getNscale(Imn)
    % Seting number of decomposition levels
    [N,M] = size(Imn);
    lmin = min(N,M);
    if lmin>550
        if lmin > 1200
            S = 3;
        else
            S=2;
        end
    else
        S=1;
    end
end



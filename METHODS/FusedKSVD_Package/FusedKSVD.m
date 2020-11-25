function [Ifd, Iss , Ims ] = FusedKSVD(Imn,NoiseSigma,DictMS,DictSS,verbose)

%  FusedKSVD: Image Denoising Algorithm through Fused K-SVD
%    
% [ Iss , Ims , Ifd ] = FusedMultiScaleDenosing(n,m,S,IterNum,DictMS,DictSS,Wname,Imn,NoiseSigma)
%   
%   This function reproduces the image denoising algorithm introduced in
%   [1], performing the fusion of a Multiscale version of the K-SVD
%   denoising algorithm and its single scale (traditional) version [2].
% 
%   Parameters:
%         Input:    - Imn: Noisy image
%                   - NoiseSigma: noise standard deviation
%                   - DictMS: Multi-scale Dictionary
%                   - DictSS: Single-scale Dictionary
% 
%         Output:   -Ifd: fused multi-scale K-SVD denoised image
%                   -Iss: single-scale K-SVD denoised image
%                   -Ims: multi-scale K-SVD denoised image
% 
% 
% [1]: J. Sulam, B. Ophir, and M. Elad, "Image Denoising Through Multi-Scale Learnt Dictionaries", 
% ICIP, Paris, October 27-30, 2014.
% 
% [2]:M. Aharon, M. Elad, and A.M. Bruckstein, "The K-SVD: An Algorithm for Designing of Overcomplete 
% Dictionaries for Sparse Representation", the IEEE Trans. On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006.

%  Note: this function requires the K-SVD and OMP box by Ron Rubinstein
%  included in the matlab path, available at: 
%  http://www.cs.technion.ac.il/~ronrubin/software.html
% 
% [1]: J. Sulam, B. Ophir, and M. Elad, "Image Denoising Through Multi-Scale Learnt Dictionaries", 
% ICIP, Paris, October 27-30, 2014.
% 
% [2]:M. Aharon, M. Elad, and A.M. Bruckstein, "The K-SVD: An Algorithm for Designing of Overcomplete 
% Dictionaries for Sparse Representation", the IEEE Trans. On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006.
%
% Jeremias Sulam
% Computer Science Department
% Technion, Haifa 3200003 Israel
% October 2014.
%
% jsulam@cs.technion.ac.il
% www.cs.technion.ac.il/~jsulam
% 
% version: 1.0
% updated: Oct.25, 2014

% 
%   checking for package
    
    if ( exist('omp2','file')==0 ||  exist('ksvd','file')==0 )
        error('ksvd and omp packages needed.')
    end

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
    
    % Wavelet name:
    Wname = 'dmey';
    
    % K-SVD iterations
    IterNum=20;
   
    [Iss, Ims ]=MScaleKSVDdenoising(S,IterNum,DictMS,DictSS,Wname,Imn,NoiseSigma,verbose);              

    % Fusing Image
    m=0.85/40;
    alfa=m*(NoiseSigma-10);
    
    if verbose == 1
        disp('Fusing...')
    end
    n = size(DictSS,1);

    Ifd = JointOMP(DictSS,Ims,Iss,alfa,sqrt(NoiseSigma^2*n)*0.09,sqrt(n));
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Iss, Ims ] = MScaleKSVDdenoising(S,IterNum,DictMS,DictSS,Wname,image,NoiseSigma,verbose)
% Multi scale and single scale K-SVD dict adaptive image denoising by K-SVD
% 
%       Jeremias Sulam
% Computer Science Dpt. - Technion
%       October 2014

TS = dwtmode('status','nodisp');
if ~strcmp(TS,'sym')
    dwtmode('sym','nodisp')
end

[n,m] = size(DictSS);

E = sqrt(NoiseSigma^2*n);   % coding threshold

[C,L]=wavedec2(image,S,Wname);

Crec=[];  

%  Regular K-SVD Denoising  %
if verbose==1
    disp('Denoising Regular K-SVD...')
end

[Iss]=CleanKSVD(image,n,m,E*getConstant(n),IterNum,DictSS,20/NoiseSigma);

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
    disp('Denoising Multi-scale K-SVD...')
end

for s=S:-1:1
    
    if s==S
        A = appcoef2(C,L,Wname,s);
        [Ac,~,~] = CleanKSVD(A,n,m,E*getConstant(n*4*s^2,0.95),IterNum,DictMs.A,20/NoiseSigma);
        Crec=[Crec,Ac(:)'];        
    end
    
    [H,V,D] = detcoef2('all',C,L,s);

    [Vc,~,~] = CleanKSVD(V,n,m,E*getConstant(n*4*s^20,0.999),IterNum,DictMs.V,0);
       
    [Hc,~,~] = CleanKSVD(H,n,m,E*getConstant(n*4*s^2,0.999),IterNum,DictMs.H,0);
    
    [Dc,~,~] = CleanKSVD(D,n,m,E*getConstant(n*4*s^2,0.999),IterNum,DictMs.D,0);
    
    Crec=[Crec,Hc(:)',Vc(:)',Dc(:)'];
        
end

Ims=waverec2(Crec,L,Wname);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageFused = JointOMP(Dict,Im1,Im2,alfa,thres,p)

% Joint OMP for image fusion with sparse coding
% 
% Parameters:
%           Dict: Dictionary
%           Im1, Im2: Images to merge
%           alfa: weighting paramter
%           threshold: coding threshold for sparse coding
%           p: patch size
% 
% OutPut:   ImageFused
% 
% Jeremias Sulam
% 01/14

[N,M] = size(Im1);
D=[Dict*sqrt(1+alfa);Dict*sqrt(1-alfa)];
D=NormDict(D);
G=D'*D;
Y=[im2col(Im1,[p,p])*sqrt(1+alfa) ; im2col(Im2,[p,p])*sqrt(1-alfa) ];

% %    subdividing for potential parallelization or reduced memory
% L = size(Y,2);
% Nparts = 4;
% X = sparse(size(D,2),L);  
% for p=1:Nparts
%     if p<Nparts
%         Yp = Y(:,1+(p-1)*floor(Y/Nparts):p*floor(Y/Nparts));
%         X(:,1+(p-1)*floor(Y/Nparts):p*floor(Y/Nparts)) = omp2(D'*Yp,sum(Yp.*Yp),G,thres);
%     else
%         Yp = Y(:,1+(p-1)*floor(Y/Nparts):end);
%         X(:,1+(p-1)*floor(Y/Nparts):end) = omp2(D'*Yp,sum(Yp.*Yp),G,thres);
%     end
% end

X = omp2(D'*Y,sum(Y.*Y),G,thres);
ImageFused = RecoverImage(N,M,Dict,X)/sqrt(2);


end

% ========================================================
% ========================================================
% ========================================================

function [yout]=RecoverImage(N,M,D,CoefMatrix)


n=sqrt(size(D,1)); 
yout=zeros(N,M); 
Weight=zeros(N,M);
i=1; j=1;
for k=1:1:(N-n+1)*(M-n+1),
    patch=reshape(D*CoefMatrix(:,k),[n,n]); 
    yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
    Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+1; 
    else
        i=1; j=j+1; 
    end;
end;
yout=(yout)./(Weight); 

end

% ========================================================
% ========================================================
% ========================================================

function D = NormDict (A)
    D=zeros(size(A));
    for i=1:size(A,2)
        if norm(A(:,i))~=0,
            D(:,i)=A(:,i)/norm(A(:,i));
        end
    end
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




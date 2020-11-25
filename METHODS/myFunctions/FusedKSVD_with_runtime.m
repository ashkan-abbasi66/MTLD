function [Ifd, Iss , Ims, runtime ] = FusedKSVD_with_runtime(Imn,NoiseSigma,DictMS,DictSS,verbose,weaken,Nscale)

% The original file is "FusedKSVD.m"
% 
% We observe that when noise level is relatively strong (such as \sigma =
% 50), the single-scale result that is obtained originally by the
% "CleanKSVD" method may have unexpected behavior. It results some dark
% spots over the denoised image. To solve this problem, we use the original
% ksvd denoising package. However, the parameters are set as close as
% possible to "CleanKSVD". This solves the mentioned problem and result in
% much faster method which have even better performance!
% 
% 
% I also added a weaken factor that can be used to shrink the threshold
% used in sparse coding stage. 
% 
% 
% Ashkan.
% 

%   checking for package
    
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
   
    [Iss, Ims , runtime] = ...
        MScaleKSVDdenoising_with_runtime(S,IterNum,DictMS,DictSS,Wname,Imn,NoiseSigma,verbose,weaken);              
    
    
    % Fusing Image
    m=0.85/40;
    alfa=m*(NoiseSigma-10);
    
    if verbose == 1
        disp('    Fusing...')
    end
    n = size(DictSS,1);
    
    tic

    Ifd = JointOMP(DictSS,Ims,Iss,alfa,sqrt(NoiseSigma^2*n)*0.09,sqrt(n));
    
    runtime.Ifd = toc;
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



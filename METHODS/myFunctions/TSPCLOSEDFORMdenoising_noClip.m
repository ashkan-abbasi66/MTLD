function [IMU,paramsout]= TSPCLOSEDFORMdenoising_noClip(I1,paramsin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See "TSPCLOSEDFORMdenoising.m"
% This function does not clip values outside the range [0,255].
%  Compared to its original version, the output of this function
%  is slightly inferior.
% Moreover, since it does not compute PSNR, it does not need the original 
%  image as input.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initializing algorithm parameters
[aa,bb]=size(I1);
sig = paramsin.sig;
iterx = paramsin.iterx;
n = paramsin.n;
N = paramsin.N;
C1 = paramsin.C;
la = paramsin.tau;
T0=  paramsin.s;
maxsp = paramsin.maxsparsity;
numiterr = paramsin.M;
if(paramsin.method==0)
    lambda0=paramsin.lambda0;
end
W =  paramsin.W;
r = paramsin.r;

threshold=C1*sig*(sqrt(n)); %\ell_2 error threshold (maximum allowed norm of the difference between a noisy patch and its denoised version) per patch

%Initial steps

%Extract image patches
[TE,idx] = my_im2col(I1,[sqrt(n),sqrt(n)],r); br=mean(TE);
TE=TE - (ones(n,1)*br); %subtract means of patches
[rows,cols] = ind2sub(size(I1)-sqrt(n)+1,idx);

N4=size(TE,2); %Total number of overlapping image patches

%Check if input training size exceeds total number of image patches
if(N4>N)
    N3=N;
else
    N3=N4;
end


de=randperm(N4);

YH=TE(:,de(1:N3));  %Use a random subset of patches in the transform learning step
STY =(ones(1,N3))*(T0);  %Vector of initial sparsity levels


%Begin iterations of the two-step denoising algorithm
for ppp=1:iterx
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l2=(lambda0)*((norm(YH,'fro'))^2); l3=l2;
    %Transform Learning Step
    if(paramsin.method==0)
        [W]= TLclosedformmethod(W,YH,numiterr,l2,l3,STY);  %Transform Learning with log-determinant+Frobenius norm regularizer
    else
        [W]= TLORTHO(W,YH,numiterr,STY);  %Orthonormal Transform Learning
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    de=randperm(N4);
    
    if(ppp==iterx)
        X1=W*TE;   %In the last iteration of the two-step denoising algorithm, apply the transform to all patches
        kT=zeros(1,size(TE,2));
    else
        YH=TE(:,de(1:N3));  %In all but the last iteration of the two-step denoising algorithm, select a random subset of patches (for use in training in the next iteration)
        X1=W*(YH);
        kT=zeros(1,size(YH,2));
    end
    
    [~,ind]=sort(abs(X1),'descend');
    er=n*(0:(size(X1,2)-1));ind=ind + (er'*ones(1,n))';
    G=(pinv([(sqrt(la)*eye(n));W])); Ga=G(:,1:n);Gb=G(:,n+1:2*n);
    
    %Variable Sparsity Update Step
    if(ppp==iterx)
        Gz=Ga*((sqrt(la))*TE);
        q=Gz;   %In the last iteration of the two-step denoising algorithm, q stores the denoised patches.
        ZS2=sqrt(sum((Gz-TE).^2));
        kT=kT+(ZS2<=threshold);  %Checks if error threshold is satisfied at zero sparsity for any of the patches
        STY=zeros(1,size(TE,2));X=zeros(n,size(TE,2)); %STY is a vector of sparsity levels and X is the corresponding sparse code matrix
    else
        Gz=Ga*((sqrt(la))*YH);
        ZS2=sqrt(sum((Gz-YH).^2));
        kT=kT+(ZS2<=threshold);  %Checks if error threshold is satisfied at zero sparsity for any of the training patches
        STY=zeros(1,size(YH,2));X=zeros(n,size(YH,2)); %STY is a vector of sparsity levels and X is the corresponding sparse code matrix
    end
    
    %Incrementing sparsity by 1 at a time, until error threshold is satisfied for all patches.
    for k=1:maxsp
        indi=find(kT==0); %Find indices of patches for which the error threshold has not yet been satisfied
        if(isempty(indi))
            break;
        end
        
        X(ind(k,indi))=X1(ind(k,indi));  %Update sparse codes to the current sparsity level in the loop.
        if(ppp==iterx)
            q(:,indi)= Gz(:,indi) + Gb*(X(:,indi));  %Update denoised patches in the last iteration of the two-step denoising algorithm
            ZS2=sqrt(sum((q(:,indi) - TE(:,indi)).^2)); kT(indi)=kT(indi)+(ZS2<=threshold);  %Check if error threshold is satisfied at sparsity k for any patches
        else
            ZS2=sqrt(sum((Gz(:,indi) + Gb*(X(:,indi)) - YH(:,indi)).^2)); kT(indi)=kT(indi)+(ZS2<=threshold); %Check if error threshold is satisfied at sparsity k for any training patches
        end
        STY(indi)=k;  %Update the sparsity levels of patches
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Averaging together the denoised patches at their respective locations in the 2D image
IMout=zeros(aa,bb);Weight=zeros(aa,bb);
bbb=sqrt(n);
for jj = 1:10000:size(TE,2)
    jumpSize = min(jj+10000-1,size(TE,2));
    ZZ= q(:,jj:jumpSize) + (ones(size(TE,1),1) * br(jj:jumpSize));
%     inx=(ZZ<0);ing= ZZ>255; ZZ(inx)=0;ZZ(ing)=255;
    for ii  = jj:jumpSize
        col = cols(ii); row = rows(ii);
        block =reshape(ZZ(:,ii-jj+1),[bbb,bbb]);
        IMout(row:row+bbb-1,col:col+bbb-1)=IMout(row:row+bbb-1,col:col+bbb-1)+block;
        Weight(row:row+bbb-1,col:col+bbb-1)=Weight(row:row+bbb-1,col:col+bbb-1)+ones(bbb);
    end
end
IMU=(IMout)./(Weight);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ing= IMU<0; ing2= IMU>255;
% IMU(ing)=0;IMU(ing2)=255;  %Limit denoised image pixel intensities to the range [0, 255].


%Output the learnt transform and the PSNR of the denoised image.
paramsout.transform=W;
% paramsout.PSNR=20*log10((sqrt(aa*bb))*255/(norm(double(IMU)-double(I7),'fro')));

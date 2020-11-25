function [ PsnrDb,MSE ] = PSNR(I,Itest,varargin)
%   [ PsnrDb,Pnsr ] = PSNR(I,Iest,R)
%   
%         input: I: true signal
%                Iest: estimated signal
%                R: max. signal value. if not specifyed, R=255;
%         output: PsnrDb: Peak signal/noise Ratio in Db
%                 Pnsr: signal/ratio in %
%                 
% Jeremias Sulam
% 05/13;

if nargin==2,
    R=255;
end
if nargin==3,
    v=cell2mat(varargin);
    R=v(1);
end

if size(I,1)~=size(Itest,1)||size(I,2)~=size(Itest,2)
    disp('Error.The data to compare must have the same size. Sorry!')
    return;
end

[n m]=size(I);

MSE=sum(sum((I-Itest).^2))/(m*n);

PsnrDb=10*log10(R^2/MSE);

end


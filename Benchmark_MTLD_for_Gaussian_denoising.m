%{


% ============
% USAGE
% ============

Only Run This Script

Set "for i = 1" to "for i = 1:N_folders", if you want to run the method 
over the whole dataset. 

If you want to calculate the metrics (PSNR, and SSIM) and 
save the results as separate image files,
after running this script, run the following piece of code.


config = {
    'Imn','','_%.3d_0_Imn'
    'Iss_tl','runtime.Iss_tl','_%.3d_4_Iss_tl';
    'Ims_tl','runtime.Ims_tl.overall','_%.3d_5_Ims_tl';   
    'Imx_dwt','runtime.Imx_dwt','_%.3d_6_Imx_dwt';
    'Imx_iuwt','runtime.Imx_iuwt','_%.3d_8_Imx_iuwt';
    };
	
args.data_path = './RESULTS/mtld_dt12';
args.save_dir = './RESULTS/mtld_dt12_imgs'; % or ''

metrics = {};

for ii = 1:size(config,1)
    args.image_variable_name = config{ii,1}; %%%%%%%%%%%%%
    args.runtime_variable_name = config{ii,2}; %%%%%%%%%%%%%
    args.postfix = config{ii,3};

    [table_psnr,table_ssim,table_time] = EVAL_Table_2(args);
    
    metrics{ii,1} = table_psnr;
    metrics{ii,2} = table_ssim;
    metrics{ii,3} = table_time;  
end


%}


%% paths

clear
clc

addpath(genpath('./METHODS'))

data_folder = 'dt12'

data_path = fullfile('./DATASETS/',data_folder); %%%%%%%%%%%%%

save_dir = fullfile('./RESULTS/',['mtld_' data_folder]); %%%%%%%%%%%%% 


%% Generate folder name for each noise level

input_sigmas = [15,25,50]; %%%%%%%%%%%%%

N_sigmas = length(input_sigmas);




%% Reading all noisy and clean images in the dataset

listing = dir(data_path);
N_folders = length(listing) - 2;

for i = 1 %  %%%%%%%%%%%%% 

    image_name = listing(i+2).name;
    image_path = fullfile(data_path,image_name);
    
    [folder_path,image_name,ext] = fileparts(image_path);
    
    msg = '\n\n====> folder %s (#%d/%d)\n';
    fprintf(msg,image_name,i,N_folders); 
    
    
    for j = 1:N_sigmas
        
        sig = input_sigmas(j);
        
        fprintf('--> noise %d (#%d/%d)\n',sig,j,N_sigmas);
        
        subfolder_name = sprintf('sigma%0.3d',sig);
        output_dir = fullfile(save_dir,image_name,subfolder_name);
        if ~exist(output_dir,'dir')
            mkdir(output_dir)
        end
        outputs_path = fullfile(output_dir,'outputs.mat');
        
        
        imn_path = fullfile(image_path,subfolder_name,'I1.mat');
        im_path = fullfile(image_path,subfolder_name,'I7.mat');
        
        load(imn_path);
        load(im_path);
        
        imn = I1(:,:,1); % noisy image
        im = I7(:,:,1); % ground-truth (clean) image
        
        
        
%% Performing denoising algorithms
        
        %=================
        % TL denoising (Iss_tl) % 
        %================= 
        n = 121; % patch size
        [Iss_tl,~,rt] = TLD(imn,im,sig,n);
        runtime.Iss_tl = rt;

        
        
        
        %=================
        % MS TL denoising (out_tl_ms.Ims_tl) %
        %=================
        wname = 'dmey';
        Ns = getNs(imn);
%         Ns = 1;
        n = 121;
        vernose = 1;
        [out_ms_tl,rt] = TLD_M(imn,sig,n,Ns,wname,vernose);
%         [out_ms_tl,rt] = TLD_M(imn,sig,n,[],wname,vernose);
        runtime.Ims_tl = rt;
        Ims_tl = out_ms_tl.Ims_tl;
        ANh = out_ms_tl.ANh;

        
        
        %=================
        % Mixed TL denoising methods (Imx_dwt, Imx_swt, Imx_iuwt) %
        %=================
        mx_wname = wname; % wavelet name for band mixing
%         mx_Ns = Ns; % number of decomposition levels for band mixing
        
       
        n = 121;
        [Imx_dwt,runtime_fmm] = ...
            TLD_FMM(imn,sig,n,Ns,Iss_tl,mx_wname); 
        runtime.Imx_dwt = runtime_fmm + runtime.Iss_tl;
        
              
        tic
        Imx_iuwt = fuse_bands_IUWT(Iss_tl,Ims_tl,Ns);
        runtime.Imx_iuwt = toc + runtime.Iss_tl + runtime.Ims_tl.overall;
        
%         imshow_eval_2(im,Iss_tl,'1.04')
%         imshow_eval_2(im,Ims_tl,'1.04 - Ims')
%         imshow_eval_2(im,Imx_dwt,'1.04 - dwt')
%         imshow_eval_2(im,Imx_iuwt,'1.04 - iuwt')

        [MM,NN] = size(im);
        ps_Iss_tl = comp_psnr_2(im,Iss_tl);
        ps_Ims_tl = comp_psnr_2(im,Ims_tl(1:MM,1:NN));
        ps_Imx_dwt = comp_psnr_2(im,Imx_dwt(1:MM,1:NN));
        ps_Imx_iuwt = comp_psnr_2(im,Imx_iuwt(1:MM,1:NN));
            
        fprintf('%9s%9s%9s%9s\n','Iss_tl','Ims_tl','Imx_dwt',...
            'Imx_iuwt')
        fprintf('%9.2f%9.2f%9.2f%9.2f\n',...
            ps_Iss_tl,ps_Ims_tl,ps_Imx_dwt,ps_Imx_iuwt)
        
        
%% Save results           
        
        %=================
        % Save MAT files %
        % =================
        
        save(outputs_path,...
            'im_path','imn_path',...            % clean and noisy image paths
            'sig',...                           % Gaussian noise level
            'Iss_tl',...
            'Ims_tl',...
            'Imx_dwt',... %Imx_swt',...
            'Imx_iuwt',...
            'runtime');  
        
        % 2:05
    end
    
end

disp('FINISHED')
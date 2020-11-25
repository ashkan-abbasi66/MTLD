%{

% ============
% USAGE
% ============

Only Run This Script

If you want to calculate the metrics (PSNR, and SSIM) and 
save the results as separate image files,
after running this script, run the following piece of code.

config = {
    'Iss','runtime.Iss','_%.3d_1_Iss';
    'Ims','runtime.Ims','_%.3d_2_Ims';
    'Ifd','runtime.Ifd','_%.3d_3_Ifd'
    };
args.data_path = './RESULTS/FusedKSVD_dt12';
args.save_dir = './RESULTS/FusedKSVD_dt12_imgs'; % or ''

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


%% Paths and Configs

clear
% close all
% clc

addpath(genpath('./METHODS'))


data_folder = 'dt12'

data_path = fullfile('./DATASETS/',data_folder); %%%%%%%%%%%%%

save_dir = fullfile('./RESULTS/',['FusedKSVD_' data_folder]); %%%%%%%%%%%%% 


%% Generate folder name for each noise level


all_sigmas = [15,25,50];

N_sigmas = length(all_sigmas);
sigma_level = cell(N_sigmas,1);
for j = 1:N_sigmas
    sigma_level{j} = sprintf('sigma%.3d',all_sigmas(j));
end



%% Reading all noisy and clean images in the dataset

listing = dir(data_path);
N_folders = length(listing) - 2;

for i = 1:N_folders

    image_folder = listing(i+2).name
    folder_path = fullfile(data_path,image_folder);
    
    msg = '\n\n====> folder %s (#%d/%d)\n';
    fprintf(msg,image_folder,i,N_folders); 
    
    
    for j = 1:N_sigmas
        
        %--------------------
        % sig = all_sigmas(j);
        str = sigma_level{j};
        str_ = str(strfind(str,'a')+1:end);
        sig = str2double(str_);
        %--------------------
        
        fprintf('--> noise %d (#%d/%d)\n',sig,j,N_sigmas);
        
        output_dir = fullfile(save_dir,image_folder,sigma_level{j});
        if ~exist(output_dir,'dir')
            mkdir(output_dir)
        end
        outputs_path = fullfile(output_dir,'outputs.mat');
        
        
        imn_path = fullfile(folder_path,sigma_level{j},'I1.mat');
        im_path = fullfile(folder_path,sigma_level{j},'I7.mat');
        
        load(imn_path);
        load(im_path);
        
        imn = I1; % noisy image
        im = I7; % ground-truth (clean) image
        
        
        
%% Performing denoising algorithms
        
        %=================
        % KSVD denoising (Iss) %
        % MS KSVD (Ims)        %
        % FusedKSVD (Ifd)      %
        %=================
        
        disp('~~ KSVD denoising methods ~~')
        
        load DictMS
        load DictSS

%         [Ifd, Iss , Ims]=FusedKSVD(imn,sig,DictMS,DictSS,1);
        [Ifd, Iss , Ims, rt] = FusedKSVD_with_runtime(imn,sig,DictMS,DictSS,1);
        
        runtime.Iss = rt.Iss;
        runtime.Ims = rt.Ims;
        runtime.Ifd = rt.Iss + rt.Ims + rt.Ifd;
        runtime.approx_denoising_time = rt.approx_denoising_time;
        
        
        
%         imshow_eval_2(im,Iss,'Iss')
%         imshow_eval_2(im,Ims,'Ims')
%         imshow_eval_2(im,Ifd,'Ifd')
        [MM,NN] = size(im);
        ps_Iss = comp_psnr_2(im,Iss);
        ps_Ims = comp_psnr_2(im,Ims(1:MM,1:NN));
        ps_Ifd = comp_psnr_2(im,Ifd(1:MM,1:NN));
            
        fprintf('%9s%9s%9s\n','Iss','Ims','Ifd')
        fprintf('%9.2f%9.2f%9.2f\n',...
            ps_Iss,ps_Ims,ps_Ifd)
        

%% Save results           
        
        %=================
        % Save MAT files %
        % =================
        
        save(outputs_path,...
            'im_path','imn_path',...            % clean and noisy image paths
            'sig',...                           % Gaussian noise level
            'Iss',...
            'Ims',...
            'Ifd',...
            'runtime');  

    end
    
end

disp('FINISHED')
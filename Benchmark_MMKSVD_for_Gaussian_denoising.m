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
    'Iss','runtime.Iss','_%.3d_1_Iss';
    'Ims','runtime.Ims','_%.3d_2_Ims';
    'Imx_dwt','runtime.Imx_dwt','_%.3d_3_Imx_dwt';
    'Imx_iuwt','runtime.Imx_iuwt','_%.3d_3_Imx_iuwt';
    };
args.data_path = './RESULTS/mmksvd_dt12';
args.save_dir = './RESULTS/mmksvd_dt12_imgs'; % or ''

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
clc

addpath(genpath('./METHODS'))


data_folder = 'dt12'; %%%%%%%%%%%%%

data_path = fullfile('./DATASETS/',data_folder); %%%%%%%%%%%%%

save_dir = fullfile('./RESULTS/',['test____' data_folder ]); %%%%%%%%%%%%% 


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

for i = 1 %%%%%%%%%%%%%

    image_folder = listing(i+2).name;
    folder_path = fullfile(data_path,image_folder);
    
    msg = '********* folder %s (#%d/%d) *********\n\n';
    fprintf(msg,image_folder,i,N_folders); 
    
    
    for j = 1:N_sigmas
        
        %--------------------
        % sig = all_sigmas(j);
        str = sigma_level{j};
        str_ = str(strfind(str,'a')+1:end);
        sig = str2double(str_);
        %--------------------
        
        fprintf('--------> noise %d (#%d/%d) \n\n',sig,j,N_sigmas);
        
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
        
      
        disp('============ KSVD denoising algorithms ============')
        
        load DictMS
        load DictSS

%         [Ifd, Iss , Ims]=FusedKSVD(imn,sig,DictMS,DictSS,1);
%         [Ifd, Iss , Ims, rt] = FusedKSVD_with_runtime(imn,sig,DictMS,DictSS,1);
        verbose = 1;
        [Iss , Ims, Imx_dwt, Imx_iuwt, rt ] = ...
            MMKSVD_with_runtime(imn,sig,DictMS,DictSS,verbose);
        
        runtime.Iss = rt.Iss;
        runtime.Ims = rt.Ims;
        runtime.Imx_dwt = rt.Imx_dwt;
        runtime.Imx_iuwt = rt.Imx_iuwt;
        
        
        
%         imshow_eval_2(im,Iss,'Iss')
%         imshow_eval_2(im,Ims,'Ims')
%         imshow_eval_2(im,Imx_dwt,'Imx-dwt')
%         imshow_eval_2(im,Imx_iuwt,'Imx-iuwt')
%% Save results           
        
        %=================
        % Save MAT files %
        % =================
        
        save(outputs_path,...
            'im_path','imn_path',...            % clean and noisy image paths
            'sig',...                           % Gaussian noise level
            'Iss',...
            'Ims','Imx_dwt','Imx_iuwt',...
            'runtime');  

    end
    
end

disp('FINISHED')
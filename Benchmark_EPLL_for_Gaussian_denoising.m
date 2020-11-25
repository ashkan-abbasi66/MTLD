%{

% ============
% USAGE
% ============


config = {
    'im_epll','runtime','_%.3d_1_epll'}
args.data_path = './RESULTS_XLSX/EPLL_dt12';
args.save_dir = './RESULTS_XLSX/EPLL_dt12_imgs'; % or ''

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

addpath(genpath('./METHODS/EPLL_code'))


data_folder = 'dt12'; %%%%%%%%%%%%%

data_path = fullfile('./DATASETS/',data_folder); %%%%%%%%%%%%%

save_dir = fullfile('./RESULTS_XLSX/',['EPLL_' data_folder ]); %%%%%%%%%%%%% 


%% Generate folder name for each noise level

% all_sigmas = [5,10,15,20,25,30,40,50,75,100]; %%%%%%%%%%%%%
all_sigmas = [15,25,50];
% all_sigmas = [50];
N_sigmas = length(all_sigmas);
sigma_level = cell(N_sigmas,1);
for j = 1:N_sigmas
    sigma_level{j} = sprintf('sigma%.3d',all_sigmas(j));
end



%% Reading all noisy and clean images in the dataset

listing = dir(data_path);
N_folders = length(listing) - 2;

for i = 1:N_folders

    image_folder = listing(i+2).name;
    folder_path = fullfile(data_path,image_folder);
    
    msg = '\n====> folder %s (#%d/%d) \n';
    fprintf(msg,image_folder,i,N_folders); 
    
    
    for j = 1:N_sigmas
        
        %--------------------
        % sig = all_sigmas(j);
        str = sigma_level{j};
        str_ = str(strfind(str,'a')+1:end);
        sig = str2double(str_);
        %--------------------
        
        fprintf('==> noise %d (#%d/%d) \n',sig,j,N_sigmas);
        
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
              
        disp('~~~ EPLL by Zoran & Weiss ')
        
        % params
        patchSize = 8;
        noiseSD = sig/255;
        excludeList = [];

        % set up prior
        LogLFunc = [];
        load GSModel_8x8_200_2M_noDC_zeromean.mat
        prior = @(Z,patchSize,noiseSD,imsize) aprxMAPGMM(Z,patchSize,noiseSD,imsize,GS,excludeList);

        % load image
        I = im/255;

        % add noise
        noiseI = imn/255;
        
        % add 64 and 128 for high noise
        tic
        [im_epll,~,~] = EPLLhalfQuadraticSplit(noiseI,patchSize^2/noiseSD^2,patchSize,(1/noiseSD^2)*[1 4 8 16 32],1,prior,I,LogLFunc);
        runtime = toc;
        
        im = im * 255;
        imn = imn * 255;
        im_epll = im_epll * 255;
        
        

%% Save results           
        
        %=================
        % Save MAT files %
        % =================
        
        save(outputs_path,...
            'im_path','imn_path',...            % clean and noisy image paths
            'sig',...                           % Gaussian noise level
            'im_epll',...
            'runtime');  

    end
    
end

disp('FINISHED')
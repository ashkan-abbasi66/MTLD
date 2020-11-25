% This script demonstrates how to use TLD and its multiscale versions 
% for removing Poisson-Guassian noise from florescence microscopy images. 
% 
% We follow the conventions introduced in the following work and its
% associated paper:
%   https://github.com/yinhaoz/denoising-fluorescence
%
% Before run this script, please download the "test_mix" subset from FMD
% dataset. 
%
% For more information, please have a look at the following link:
%   https://github.com/ashkan-abbasi66
%
%  ========== DATE ==========
%  9/21/2020



%% Paths and Configs

clear
clc

addpath(genpath('./METHODS/'))

dataset = 'test_mix';
dataset_dir = fullfile('./DATASETS/',dataset); 

% Indices of input images for denoising
which_images = 37:48; % Set to 1:48 for selecting all images in `test_mix`
N_images = length(which_images); % Number of images


iuwt_func = @fuse_bands_IUWT_79;


result_path = './RESULTS/';
dir_name = 'test_mix_TLD98_MTLD9_37_48';
save_dir = fullfile(result_path,dir_name);


sp_1 = 0.98; % sparsity controller for TLD
sp_2 = 0.9; % sparsity controller for MTLD
Nscale = 1;


% MAT-file to store all results
out_matfile = fullfile(result_path,[dir_name '.mat']);


% Read images from here
dir_avg1 = fullfile(dataset_dir, '/raw/'); 
dir_avg2 = fullfile(dataset_dir, '/avg2/');
dir_avg4 = fullfile(dataset_dir, '/avg4/');
dir_avg8 = fullfile(dataset_dir, '/avg8/'); 
dir_avg16 = fullfile(dataset_dir, '/avg16/');
dir_gt = fullfile(dataset_dir, '/gt/');


% Save results here
save_dir1 = fullfile(save_dir,'avg1'); 
save_dir2 = fullfile(save_dir,'avg2');
save_dir4 = fullfile(save_dir,'avg4');
save_dir8 = fullfile(save_dir,'avg8');
save_dir16 = fullfile(save_dir,'avg16');

% Output images are saved with the following postfixes in their names
postfix_imn = '_1_noisy';
postfix_imout_tld = '_2_TLD';
postfix_imout_mtld = '_3_MTLD';
postfix_imout_fmmtld = '_4_FMMTLD';
postfix_imout_mmtld = '_5_MMTLD';



%% Import images and create arrays to store PSNR and SSIM values

N_methods = 5;

[avg1_array,img_names] = import_img_array_2(dir_avg1);
img_names = img_names(which_images);

avg1_PSNRs = zeros(N_images,N_methods);
avg1_SSIMs = zeros(N_images,N_methods);
avg1_TIMEs = zeros(N_images,1);

[avg2_array] = import_img_array(dir_avg2);
avg2_PSNRs = zeros(N_images,N_methods);
avg2_SSIMs = zeros(N_images,N_methods);
avg2_TIMEs = zeros(N_images,1);


[avg4_array] = import_img_array(dir_avg4);
avg4_PSNRs = zeros(N_images,N_methods);
avg4_SSIMs = zeros(N_images,N_methods);
avg4_TIMEs = zeros(N_images,1);


[avg8_array] = import_img_array(dir_avg8);
avg8_PSNRs = zeros(N_images,N_methods);
avg8_SSIMs = zeros(N_images,N_methods);
avg8_TIMEs = zeros(N_images,1);

[avg16_array] = import_img_array(dir_avg16);
avg16_PSNRs = zeros(N_images,N_methods);
avg16_SSIMs = zeros(N_images,N_methods);
avg16_TIMEs = zeros(N_images,1);

% Import ground-truth images
[img_gt_array] = import_img_array(dir_gt);



%% Create cell arrays to store all denoised images

titles = {'file name','imn',...
    'imout_tld','imout_tld_direct',...
    'imout_mtld','mtld_outputs',...
    'imout_fmmtld','imout_mmtld'};

denoised_avg1 = cell(length(img_names)+1,length(titles));
denoised_avg1(1,:) = titles; % first row
denoised_avg1(2:N_images+1,1) = img_names; % first column

denoised_avg2 = cell(length(img_names)+1,length(titles));
denoised_avg2(1,:) = titles;
denoised_avg2(2:N_images+1,1) = img_names;

denoised_avg4 = cell(length(img_names)+1,length(titles));
denoised_avg4(1,:) = titles;
denoised_avg4(2:N_images+1,1) = img_names;

denoised_avg8 = cell(length(img_names)+1,length(titles));
denoised_avg8(1,:) = titles;
denoised_avg8(2:N_images+1,1) = img_names;

denoised_avg16 = cell(length(img_names)+1,length(titles));
denoised_avg16(1,:) = titles;
denoised_avg16(2:N_images+1,1) = img_names;



for counter = 1:length(which_images)
        
    fprintf('\n----- IMAGE %d / %d -----\n', counter,N_images)

    img_name = img_names{counter};
    [pathstr, name, ext] = fileparts(img_name);
    img_name = name;

    fprintf('-- %s\n',img_name)
    disp('')
    
    img_index = which_images(counter);


    % Read 5 images with different noise realizations

    img_avg1 = avg1_array(:,:,img_index);
    img_avg2 = avg2_array(:,:,img_index);
    img_avg4 = avg4_array(:,:,img_index);
    img_avg8 = avg8_array(:,:,img_index);
    img_avg16 = avg16_array(:,:,img_index);   

    % Read ground-truth image
    img_gt = img_gt_array(:,:,img_index);% 


    
    %% Denoising ...


    %% == AVG1 ==
    
    disp('AVG1')
    
    % imout_tld
    [imout_tld, runtime_iss, imout_tld_direct] = denoise_VST_TLD(img_avg1,sp_1);
    rt_TLD = runtime_iss;
    
    % imout_mtld
    [imout_mtld, runtime_ims, mtld_outputs] = denoise_VST_MTLD(img_avg1,sp_2,Nscale);
    rt_MTLD = runtime_ims.overall;

    % ~~~ FMMTLD ~~~~
    % To obtian the FMMTLD result, you can use one of the followings:
    % The following two methods should give similar results. However, I see
    % that the second one is slightly better, especially when noise becomes
    % weak.
%     [imout_fmmtld, rt_fmmtld] = denoise_VST_FMMTLD(img_avg1,sp_2,Nscale,imout_tld);
%     rt_dwt = rt_fmmtld + rt_TLD;
    
    
    tic
    [imout_fmmtld] = change_A_band_and_inverse_Anscombe(imout_tld_direct,mtld_outputs);
    rt_dwt = toc + rt_TLD + runtime_ims.ANh;
    % ~~~~~~~~~~~~~~~~

    
    % imout_mmtld
    tic
    [imout_mmtld] = fuse_bands_IUWT_and_inverse_Anscombe(...
        imout_tld_direct,mtld_outputs,iuwt_func);
    rt_iuwt = toc + rt_TLD + runtime_ims.overall;


    % To fill "denoised_avgX" cell array, see "titles"
    denoised_avg1{img_index + 1,2} = img_avg1; % ***********
    denoised_avg1{img_index + 1,3} = imout_tld; % ***********
    denoised_avg1{img_index + 1,4} = imout_tld_direct;
    denoised_avg1{img_index + 1,5} = imout_mtld; % ***********
    denoised_avg1{img_index + 1,6} = mtld_outputs;
    denoised_avg1{img_index + 1,7} = imout_fmmtld; % ***********
    denoised_avg1{img_index + 1,8} = imout_mmtld; % ***********

    
    % Store runtimes
    
    avg1_TIMEs(img_index,1) = 0;
    avg1_TIMEs(img_index,2) = rt_TLD;
    avg1_TIMEs(img_index,3) = rt_MTLD;
    avg1_TIMEs(img_index,4) = rt_dwt;
    avg1_TIMEs(img_index,5) = rt_iuwt;
    
    
    
    %% == AVG2 ==
    
    disp('AVG2')

    % imout_tld
    [imout_tld, ~,imout_tld_direct] = denoise_VST_TLD(img_avg2,sp_1);
    rt_TLD = runtime_iss;
    
    
    % imout_mtld
    [imout_mtld, ~,mtld_outputs] = denoise_VST_MTLD(img_avg2,sp_2);
    rt_MTLD = runtime_ims.overall;
    
    
    % imout_fmmtld  
    tic
    [imout_fmmtld] = change_A_band_and_inverse_Anscombe(imout_tld_direct,mtld_outputs);
    rt_dwt = toc + rt_TLD + runtime_ims.ANh;
    
    
    % imout_mmtld
    tic
    [imout_mmtld] = fuse_bands_IUWT_and_inverse_Anscombe(...
        imout_tld_direct,mtld_outputs,iuwt_func);
    rt_iuwt = toc + rt_TLD + runtime_ims.overall;
    

    denoised_avg2{img_index + 1,2} = img_avg2; 
    denoised_avg2{img_index + 1,3} = imout_tld;
    denoised_avg2{img_index + 1,4} = imout_tld_direct;
    denoised_avg2{img_index + 1,5} = imout_mtld;
    denoised_avg2{img_index + 1,6} = mtld_outputs;
    denoised_avg2{img_index + 1,7} = imout_fmmtld;
    denoised_avg2{img_index + 1,8} = imout_mmtld;


    % Store runtimes
    
    avg2_TIMEs(img_index,1) = 0;
    avg2_TIMEs(img_index,2) = rt_TLD;
    avg2_TIMEs(img_index,3) = rt_MTLD;
    avg2_TIMEs(img_index,4) = rt_dwt;
    avg2_TIMEs(img_index,5) = rt_iuwt;
    
    
    %% == AVG4 ==

    disp('AVG4')
    
    % imout_tld
    [imout_tld, ~,imout_tld_direct] = denoise_VST_TLD(img_avg4,sp_1);
    rt_TLD = runtime_iss;
    
    
    % imout_mtld
    [imout_mtld, ~,mtld_outputs] = denoise_VST_MTLD(img_avg4,sp_2);
    rt_MTLD = runtime_ims.overall;

    
    % imout_fmmtld
    tic
    [imout_fmmtld] = change_A_band_and_inverse_Anscombe(imout_tld_direct,mtld_outputs);
    rt_dwt = toc + rt_TLD + runtime_ims.ANh;
    
    
    % imout_mmtld
    tic
    [imout_mmtld] = fuse_bands_IUWT_and_inverse_Anscombe(...
        imout_tld_direct,mtld_outputs,iuwt_func);
    rt_iuwt = toc + rt_TLD + runtime_ims.overall;


    denoised_avg4{img_index + 1,2} = img_avg4; 
    denoised_avg4{img_index + 1,3} = imout_tld;
    denoised_avg4{img_index + 1,4} = imout_tld_direct;
    denoised_avg4{img_index + 1,5} = imout_mtld;
    denoised_avg4{img_index + 1,6} = mtld_outputs;
    denoised_avg4{img_index + 1,7} = imout_fmmtld;
    denoised_avg4{img_index + 1,8} = imout_mmtld;

    
    
    % Store runtimes
    
    avg4_TIMEs(img_index,1) = 0;
    avg4_TIMEs(img_index,2) = rt_TLD;
    avg4_TIMEs(img_index,3) = rt_MTLD;
    avg4_TIMEs(img_index,4) = rt_dwt;
    avg4_TIMEs(img_index,5) = rt_iuwt;
    
    

    %% == AVG8 ==
    
    disp('AVG8')

    % imout_tld
    [imout_tld, ~,imout_tld_direct] = denoise_VST_TLD(img_avg8,sp_1);
    rt_TLD = runtime_iss;
    
    
    % imout_mtld
    [imout_mtld, ~,mtld_outputs] = denoise_VST_MTLD(img_avg8,sp_2);
    rt_MTLD = runtime_ims.overall;

    
    % imout_fmmtld
    tic
    [imout_fmmtld] = change_A_band_and_inverse_Anscombe(imout_tld_direct,mtld_outputs);
    rt_dwt = toc + rt_TLD + runtime_ims.ANh;
    
    
    % imout_mmtld
    tic
    [imout_mmtld] = fuse_bands_IUWT_and_inverse_Anscombe(...
        imout_tld_direct,mtld_outputs,iuwt_func);
    rt_iuwt = toc + rt_TLD + runtime_ims.overall;


    denoised_avg8{img_index + 1,2} = img_avg8; 
    denoised_avg8{img_index + 1,3} = imout_tld;
    denoised_avg8{img_index + 1,4} = imout_tld_direct;
    denoised_avg8{img_index + 1,5} = imout_mtld;
    denoised_avg8{img_index + 1,6} = mtld_outputs;
    denoised_avg8{img_index + 1,7} = imout_fmmtld;
    denoised_avg8{img_index + 1,8} = imout_mmtld;

    % Store runtimes
    
    avg8_TIMEs(img_index,1) = 0;
    avg8_TIMEs(img_index,2) = rt_TLD;
    avg8_TIMEs(img_index,3) = rt_MTLD;
    avg8_TIMEs(img_index,4) = rt_dwt;
    avg8_TIMEs(img_index,5) = rt_iuwt;

    
    
    %% == AVG16 ==
    
    disp('AVG16')
    
    % imout_tld
    [imout_tld, ~,imout_tld_direct] = denoise_VST_TLD(img_avg16,sp_1);
    rt_TLD = runtime_iss;
    
    
    % imout_mtld
    [imout_mtld, ~,mtld_outputs] = denoise_VST_MTLD(img_avg16,sp_2);
    rt_MTLD = runtime_ims.overall;
    

    % imout_fmmtld
    tic
    [imout_fmmtld] = change_A_band_and_inverse_Anscombe(imout_tld_direct,mtld_outputs);
    rt_dwt = toc + rt_TLD + runtime_ims.ANh;

    
    % imout_mmtld
    tic
    [imout_mmtld] = fuse_bands_IUWT_and_inverse_Anscombe(...
        imout_tld_direct,mtld_outputs,iuwt_func);
    rt_iuwt = toc + rt_TLD + runtime_ims.overall;

    
    denoised_avg16{img_index + 1,2} = img_avg16; 
    denoised_avg16{img_index + 1,3} = imout_tld;
    denoised_avg16{img_index + 1,4} = imout_tld_direct;
    denoised_avg16{img_index + 1,5} = imout_mtld;
    denoised_avg16{img_index + 1,6} = mtld_outputs;
    denoised_avg16{img_index + 1,7} = imout_fmmtld;
    denoised_avg16{img_index + 1,8} = imout_mmtld;

    % Store runtimes
    
    avg16_TIMEs(img_index,1) = 0;
    avg16_TIMEs(img_index,2) = rt_TLD;
    avg16_TIMEs(img_index,3) = rt_MTLD;
    avg16_TIMEs(img_index,4) = rt_dwt;
    avg16_TIMEs(img_index,5) = rt_iuwt;
    
    
    
    %% PSNRs
    
    disp('Compute Metrics ...')
    
    % =============
    % AVG1
    % =============

    avg1_PSNRs(img_index,1) = comp_psnr_2(img_gt,img_avg1);
    avg1_PSNRs(img_index,2) = comp_psnr_2(img_gt,denoised_avg1{img_index + 1,3}); % imout_tld
    avg1_PSNRs(img_index,3) = comp_psnr_2(img_gt,denoised_avg1{img_index + 1,5}); % imout_mtld
    avg1_PSNRs(img_index,4) = comp_psnr_2(img_gt,denoised_avg1{img_index + 1,7}); % imout_fmmtld
    avg1_PSNRs(img_index,5) = comp_psnr_2(img_gt,denoised_avg1{img_index + 1,8}); % imout_mmtld


    % =============
    % AVG2
    % =============

    avg2_PSNRs(img_index,1) = comp_psnr_2(img_gt,img_avg2);
    avg2_PSNRs(img_index,2) = comp_psnr_2(img_gt,denoised_avg2{img_index + 1,3}); % imout_tld
    avg2_PSNRs(img_index,3) = comp_psnr_2(img_gt,denoised_avg2{img_index + 1,5}); % imout_mtld
    avg2_PSNRs(img_index,4) = comp_psnr_2(img_gt,denoised_avg2{img_index + 1,7}); % imout_fmmtld
    avg2_PSNRs(img_index,5) = comp_psnr_2(img_gt,denoised_avg2{img_index + 1,8}); % imout_mmtld

    % =============
    % AVG4
    % =============

    avg4_PSNRs(img_index,1) = comp_psnr_2(img_gt,img_avg4);
    avg4_PSNRs(img_index,2) = comp_psnr_2(img_gt,denoised_avg4{img_index + 1,3}); % imout_tld
    avg4_PSNRs(img_index,3) = comp_psnr_2(img_gt,denoised_avg4{img_index + 1,5}); % imout_mtld
    avg4_PSNRs(img_index,4) = comp_psnr_2(img_gt,denoised_avg4{img_index + 1,7}); % imout_fmmtld
    avg4_PSNRs(img_index,5) = comp_psnr_2(img_gt,denoised_avg4{img_index + 1,8}); % imout_mmtld

    % =============
    % AVG8
    % =============

    avg8_PSNRs(img_index,1) = comp_psnr_2(img_gt,img_avg8);
    avg8_PSNRs(img_index,2) = comp_psnr_2(img_gt,denoised_avg8{img_index + 1,3}); % imout_tld
    avg8_PSNRs(img_index,3) = comp_psnr_2(img_gt,denoised_avg8{img_index + 1,5}); % imout_mtld
    avg8_PSNRs(img_index,4) = comp_psnr_2(img_gt,denoised_avg8{img_index + 1,7}); % imout_fmmtld
    avg8_PSNRs(img_index,5) = comp_psnr_2(img_gt,denoised_avg8{img_index + 1,8}); % imout_mmtld


    % =============
    % AVG16
    % =============

    avg16_PSNRs(img_index,1) = comp_psnr_2(img_gt,img_avg16);
    avg16_PSNRs(img_index,2) = comp_psnr_2(img_gt,denoised_avg16{img_index + 1,3}); % imout_tld
    avg16_PSNRs(img_index,3) = comp_psnr_2(img_gt,denoised_avg16{img_index + 1,5}); % imout_mtld
    avg16_PSNRs(img_index,4) = comp_psnr_2(img_gt,denoised_avg16{img_index + 1,7}); % imout_fmmtld
    avg16_PSNRs(img_index,5) = comp_psnr_2(img_gt,denoised_avg16{img_index + 1,8}); % imout_mmtld



    %% SSIMs


    avg1_SSIMs(img_index,1) = ssim(img_avg1,img_gt);
    avg1_SSIMs(img_index,2) = ssim(denoised_avg1{img_index + 1,3},img_gt); % imout_tld
    avg1_SSIMs(img_index,3) = ssim(denoised_avg1{img_index + 1,5},img_gt); % imout_mtld
    avg1_SSIMs(img_index,4) = ssim(denoised_avg1{img_index + 1,7},img_gt); % imout_fmmtld
    avg1_SSIMs(img_index,5) = ssim(denoised_avg1{img_index + 1,8},img_gt); % imout_mmtld


    avg2_SSIMs(img_index,1) = ssim(img_avg2,img_gt);
    avg2_SSIMs(img_index,2) = ssim(denoised_avg2{img_index + 1,3},img_gt); % imout_tld
    avg2_SSIMs(img_index,3) = ssim(denoised_avg2{img_index + 1,5},img_gt); % imout_mtld
    avg2_SSIMs(img_index,4) = ssim(denoised_avg2{img_index + 1,7},img_gt); % imout_fmmtld
    avg2_SSIMs(img_index,5) = ssim(denoised_avg2{img_index + 1,8},img_gt); % imout_mmtld


    avg4_SSIMs(img_index,1) = ssim(img_avg4,img_gt);
    avg4_SSIMs(img_index,2) = ssim(denoised_avg4{img_index + 1,3},img_gt); % imout_tld
    avg4_SSIMs(img_index,3) = ssim(denoised_avg4{img_index + 1,5},img_gt); % imout_mtld
    avg4_SSIMs(img_index,4) = ssim(denoised_avg4{img_index + 1,7},img_gt); % imout_fmmtld
    avg4_SSIMs(img_index,5) = ssim(denoised_avg4{img_index + 1,8},img_gt); % imout_mmtld


    avg8_SSIMs(img_index,1) = ssim(img_avg8,img_gt);
    avg8_SSIMs(img_index,2) = ssim(denoised_avg8{img_index + 1,3},img_gt); % imout_tld
    avg8_SSIMs(img_index,3) = ssim(denoised_avg8{img_index + 1,5},img_gt); % imout_mtld
    avg8_SSIMs(img_index,4) = ssim(denoised_avg8{img_index + 1,7},img_gt); % imout_fmmtld
    avg8_SSIMs(img_index,5) = ssim(denoised_avg8{img_index + 1,8},img_gt); % imout_mmtld


    avg16_SSIMs(img_index,1) = ssim(img_avg16,img_gt);
    avg16_SSIMs(img_index,2) = ssim(denoised_avg16{img_index + 1,3},img_gt); % imout_tld
    avg16_SSIMs(img_index,3) = ssim(denoised_avg16{img_index + 1,5},img_gt); % imout_mtld
    avg16_SSIMs(img_index,4) = ssim(denoised_avg16{img_index + 1,7},img_gt); % imout_fmmtld
    avg16_SSIMs(img_index,5) = ssim(denoised_avg16{img_index + 1,8},img_gt); % imout_mmtld



    %% Save Images
    
        if ~isempty(save_dir)
            
            disp('Save Images ...')

            if ~exist(save_dir1,'dir')
                mkdir(save_dir1)
                mkdir(save_dir2)
                mkdir(save_dir4)
                mkdir(save_dir8)
                mkdir(save_dir16)
            end


            %===================
            % Save output image
            %===================
            ext_ = '.tif';


            %
            % AVG1
            %
            name_str = [img_name postfix_imn ext_];
            imwrite(denoised_avg1{img_index + 1,2},fullfile(save_dir1,name_str))

            name_str = [img_name postfix_imout_tld ext_];
            imwrite(denoised_avg1{img_index + 1,3},fullfile(save_dir1,name_str))

            name_str = [img_name postfix_imout_mtld ext_];
            imwrite(denoised_avg1{img_index + 1,5},fullfile(save_dir1,name_str))

            name_str = [img_name postfix_imout_fmmtld ext_];
            imwrite(denoised_avg1{img_index + 1,7},fullfile(save_dir1,name_str))

            name_str = [img_name postfix_imout_mmtld ext_];
            imwrite(denoised_avg1{img_index + 1,8},fullfile(save_dir1,name_str))



            %
            % AVG2
            %
            name_str = [img_name postfix_imn ext_];
            imwrite(denoised_avg2{img_index + 1,2},fullfile(save_dir2,name_str))

            name_str = [img_name postfix_imout_tld ext_];
            imwrite(denoised_avg2{img_index + 1,3},fullfile(save_dir2,name_str))

            name_str = [img_name postfix_imout_mtld ext_];
            imwrite(denoised_avg2{img_index + 1,5},fullfile(save_dir2,name_str))

            name_str = [img_name postfix_imout_fmmtld ext_];
            imwrite(denoised_avg2{img_index + 1,7},fullfile(save_dir2,name_str))

            name_str = [img_name postfix_imout_mmtld ext_];
            imwrite(denoised_avg2{img_index + 1,8},fullfile(save_dir2,name_str))



            %
            % AVG4
            %
            name_str = [img_name postfix_imn ext_];
            imwrite(denoised_avg4{img_index + 1,2},fullfile(save_dir4,name_str))

            name_str = [img_name postfix_imout_tld ext_];
            imwrite(denoised_avg4{img_index + 1,3},fullfile(save_dir4,name_str))

            name_str = [img_name postfix_imout_mtld ext_];
            imwrite(denoised_avg4{img_index + 1,5},fullfile(save_dir4,name_str))

            name_str = [img_name postfix_imout_fmmtld ext_];
            imwrite(denoised_avg4{img_index + 1,7},fullfile(save_dir4,name_str))

            name_str = [img_name postfix_imout_mmtld ext_];
            imwrite(denoised_avg4{img_index + 1,8},fullfile(save_dir4,name_str))



            %
            % AVG8
            %
            name_str = [img_name postfix_imn ext_];
            imwrite(denoised_avg8{img_index + 1,2},fullfile(save_dir8,name_str))

            name_str = [img_name postfix_imout_tld ext_];
            imwrite(denoised_avg8{img_index + 1,3},fullfile(save_dir8,name_str))

            name_str = [img_name postfix_imout_mtld ext_];
            imwrite(denoised_avg8{img_index + 1,5},fullfile(save_dir8,name_str))

            name_str = [img_name postfix_imout_fmmtld ext_];
            imwrite(denoised_avg8{img_index + 1,7},fullfile(save_dir8,name_str))

            name_str = [img_name postfix_imout_mmtld ext_];
            imwrite(denoised_avg8{img_index + 1,8},fullfile(save_dir8,name_str))



            %
            % AVG16
            %
            name_str = [img_name postfix_imn ext_];
            imwrite(denoised_avg16{img_index + 1,2},fullfile(save_dir16,name_str))

            name_str = [img_name postfix_imout_tld ext_];
            imwrite(denoised_avg16{img_index + 1,3},fullfile(save_dir16,name_str))

            name_str = [img_name postfix_imout_mtld ext_];
            imwrite(denoised_avg16{img_index + 1,5},fullfile(save_dir16,name_str))

            name_str = [img_name postfix_imout_fmmtld ext_];
            imwrite(denoised_avg16{img_index + 1,7},fullfile(save_dir16,name_str))

            name_str = [img_name postfix_imout_mmtld ext_];
            imwrite(denoised_avg16{img_index + 1,8},fullfile(save_dir16,name_str))


            %{ 
            % Save ground-truth image

            to_save_fpath = fullfile(save_dir,[img_name ext_]);
            if ~exist(to_save_fpath,'file')
                imwrite(im/maxI,to_save_fpath)
            end
            %}
        end

end

%% Save a .MAT file to store all results

disp('Save MAT-file ...')

save(out_matfile,...
'img_names',...
'denoised_avg1','denoised_avg2',...
'denoised_avg4','denoised_avg8',...
'denoised_avg16',...
'avg1_PSNRs','avg1_SSIMs',...
'avg2_PSNRs','avg2_SSIMs',...
'avg4_PSNRs','avg4_SSIMs',...
'avg8_PSNRs','avg8_SSIMs',...
'avg16_PSNRs','avg16_SSIMs',...
'avg1_TIMEs','avg2_TIMEs','avg4_TIMEs',...
'avg8_TIMEs','avg16_TIMEs');


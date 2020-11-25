%{


9/22/2020
%}



%%
clc
clear

addpath(genpath('./METHODS/'))




%% INPUT: MIXED DATASET 

dataset = 'test_mix';
dataset_dir = fullfile('./DATASETS/',dataset); 


% Indices of input images for denoising
which_images = 13:16%1:48; % Set to 1:48 for selecting all images in `test_mix`
N_images = length(which_images); % Number of images


result_path = './RESULTS_XLSX/test_mix/';
dir_name = 'test_mix_FusedKSVD_13_16_final_newer';
save_dir = fullfile(result_path,dir_name);


mat_fpath = fullfile(result_path,[dir_name '.mat']);

% log_file_path = fullfile(result_path,[dir_name '.log']);
% diary(log_file_path)


dataset_path = 'E:/DTSET/FMDD/test_mix/';


dir_avg1 = [dataset_path, '/raw/']; 
dir_avg2 = [dataset_path, '/avg2/'];
dir_avg4 = [dataset_path, '/avg4/'];
dir_avg8 = [dataset_path, '/avg8/']; 
dir_avg16 = [dataset_path, '/avg16/'];
dir_gt = [dataset_path, '/gt/'];


save_dir1 = fullfile(save_dir,'avg1'); 
save_dir2 = fullfile(save_dir,'avg2');
save_dir4 = fullfile(save_dir,'avg4');
save_dir8 = fullfile(save_dir,'avg8');
save_dir16 = fullfile(save_dir,'avg16');

postfix_Iss = '_1_Iss';
postfix_Ims = '_2_Ims';
postfix_Ifd = '_3_Ifd';


%% Import images and create arrays to store PSNR and SSIM values

N_methods = 4; % including one column for the noisy image

[avg1_array,img_names] = import_img_array_2(dir_avg1);
img_names = img_names(which_images);

avg1_PSNRs = zeros(N_images,N_methods);
avg1_SSIMs = zeros(N_images,N_methods);
avg1_TIMEs = zeros(N_images,N_methods);

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
    'iss','ims',...
    'ifd'};

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



%% Loop over input images to run the denoising methods

for ii = 1:length(which_images)
    
    img_index = which_images(ii);
    
    fprintf('\n----- IMAGE %d / %d -----\n', ii,N_images)

    img_name = img_names{ii};
    [pathstr, name, ext] = fileparts(img_name);
    img_name = name;

    fprintf('----- %s\n',img_name)
    disp('')
    
    
    % Read 5 images with different noise realizations
    
    img_avg1 = avg1_array(:,:,img_index);
    img_avg2 = avg2_array(:,:,img_index);
    img_avg4 = avg4_array(:,:,img_index);
    img_avg8 = avg8_array(:,:,img_index);
    img_avg16 = avg16_array(:,:,img_index);   
    
    img_gt = img_gt_array(:,:,img_index);% for INPUT 2
    
   
    %% Perform denoising algorithms
    
    %% == AVG1 ==

    
    disp('AVG1')
    
    [Ifd_avg1, Iss_avg1 , Ims_avg1, runtime] = denoise_VST_FusedKSVD(img_avg1);
    % To fill "denoised_avgX" cell array, see "titles"
    denoised_avg1{img_index + 1,2} = img_avg1;
    denoised_avg1{img_index + 1,3} = Iss_avg1;
    denoised_avg1{img_index + 1,4} = Ims_avg1;
    denoised_avg1{img_index + 1,5} = Ifd_avg1;
    
    % Store runtimes
    
    avg1_TIMEs(img_index,1:N_methods) = [0,runtime.Iss,runtime.Ims,runtime.Ifd];
    
    
    %% == AVG2 == 
    
    disp('AVG2')
    
    [Ifd_avg2, Iss_avg2 , Ims_avg2, runtime] = denoise_VST_FusedKSVD(img_avg2);
    denoised_avg2{img_index + 1,2} = img_avg2;
    denoised_avg2{img_index + 1,3} = Iss_avg2;
    denoised_avg2{img_index + 1,4} = Ims_avg2;
    denoised_avg2{img_index + 1,5} = Ifd_avg2;
    
    % Store runtimes
    
    avg2_TIMEs(img_index,1:N_methods) = [0,runtime.Iss,runtime.Ims,runtime.Ifd];
    
    %% == AVG4 ==  
    
    
    disp('AVG4')
    
    [Ifd_avg4, Iss_avg4 , Ims_avg4, runtime] = denoise_VST_FusedKSVD(img_avg4);
    denoised_avg4{img_index + 1,2} = img_avg4;
    denoised_avg4{img_index + 1,3} = Iss_avg4;
    denoised_avg4{img_index + 1,4} = Ims_avg4;
    denoised_avg4{img_index + 1,5} = Ifd_avg4;
    
    % Store runtimes
    
    avg4_TIMEs(img_index,1:N_methods) = [0,runtime.Iss,runtime.Ims,runtime.Ifd];
    
    
    %% == AVG8 == 
    
    
    disp('AVG8')
    
    [Ifd_avg8, Iss_avg8 , Ims_avg8, runtime] = denoise_VST_FusedKSVD(img_avg8);
    denoised_avg8{img_index + 1,2} = img_avg8;
    denoised_avg8{img_index + 1,3} = Iss_avg8;
    denoised_avg8{img_index + 1,4} = Ims_avg8;
    denoised_avg8{img_index + 1,5} = Ifd_avg8;
    
    % Store runtimes
    
    avg8_TIMEs(img_index,1:N_methods) = [0,runtime.Iss,runtime.Ims,runtime.Ifd];
    
    
    %% == AVG16 == 
    % =============
    
    disp('AVG16')
    
    [Ifd_avg16, Iss_avg16 , Ims_avg16, runtime] = denoise_VST_FusedKSVD(img_avg16);
    denoised_avg16{img_index + 1,2} = img_avg16;
    denoised_avg16{img_index + 1,3} = Iss_avg16;
    denoised_avg16{img_index + 1,4} = Ims_avg16;
    denoised_avg16{img_index + 1,5} = Ifd_avg16;

    % Store runtimes
    
    avg16_TIMEs(img_index,1:N_methods) = [0,runtime.Iss,runtime.Ims,runtime.Ifd];
    
    

%% Compute PSNRs - img_gt: ground-truth (clean) image

    disp('Compute Metrics ...')
    
    maxI = 1;
    
    avg1_PSNRs(ii,1) = comp_psnr_2(img_gt,img_avg1);
    avg1_PSNRs(ii,2) = comp_psnr_2(img_gt,Iss_avg1);
    avg1_PSNRs(ii,3) = comp_psnr_2(img_gt,Ims_avg1);
    avg1_PSNRs(ii,4) = comp_psnr_2(img_gt,Ifd_avg1);
    
    
    avg2_PSNRs(ii,1) = comp_psnr_2(img_gt,img_avg2);
    avg2_PSNRs(ii,2) = comp_psnr_2(img_gt,Iss_avg2);
    avg2_PSNRs(ii,3) = comp_psnr_2(img_gt,Ims_avg2);
    avg2_PSNRs(ii,4) = comp_psnr_2(img_gt,Ifd_avg2);
    
    avg4_PSNRs(ii,1) = comp_psnr_2(img_gt,img_avg4);
    avg4_PSNRs(ii,2) = comp_psnr_2(img_gt,Iss_avg4);
    avg4_PSNRs(ii,3) = comp_psnr_2(img_gt,Ims_avg4);
    avg4_PSNRs(ii,4) = comp_psnr_2(img_gt,Ifd_avg4);
    
    avg8_PSNRs(ii,1) = comp_psnr_2(img_gt,img_avg8);
    avg8_PSNRs(ii,2) = comp_psnr_2(img_gt,Iss_avg8);
    avg8_PSNRs(ii,3) = comp_psnr_2(img_gt,Ims_avg8);
    avg8_PSNRs(ii,4) = comp_psnr_2(img_gt,Ifd_avg8);
    
    avg16_PSNRs(ii,1) = comp_psnr_2(img_gt,img_avg16);
    avg16_PSNRs(ii,2) = comp_psnr_2(img_gt,Iss_avg16);
    avg16_PSNRs(ii,3) = comp_psnr_2(img_gt,Ims_avg16);
    avg16_PSNRs(ii,4) = comp_psnr_2(img_gt,Ifd_avg16);

    
    
%% Compute SSIM    
    
    [~,avg1_SSIMs(ii,1)] = comp_psnr_2(img_avg1, img_gt);
    [~,avg1_SSIMs(ii,2)] = comp_psnr_2(Iss_avg1, img_gt);
    [~,avg1_SSIMs(ii,3)] = comp_psnr_2(Ims_avg1, img_gt);
    [~,avg1_SSIMs(ii,4)] = comp_psnr_2(Ifd_avg1, img_gt);

    [~,avg2_SSIMs(ii,1)] = comp_psnr_2(img_avg2, img_gt);
    [~,avg2_SSIMs(ii,2)] = comp_psnr_2(Iss_avg2, img_gt);
    [~,avg2_SSIMs(ii,3)] = comp_psnr_2(Ims_avg2, img_gt);
    [~,avg2_SSIMs(ii,4)] = comp_psnr_2(Ifd_avg2, img_gt);
    
    [~,avg4_SSIMs(ii,1)] = comp_psnr_2(img_avg4, img_gt);
    [~,avg4_SSIMs(ii,2)] = comp_psnr_2(Iss_avg4, img_gt);
    [~,avg4_SSIMs(ii,3)] = comp_psnr_2(Ims_avg4, img_gt);
    [~,avg4_SSIMs(ii,4)] = comp_psnr_2(Ifd_avg4, img_gt);

    [~,avg8_SSIMs(ii,1)] = comp_psnr_2(img_avg8, img_gt);
    [~,avg8_SSIMs(ii,2)] = comp_psnr_2(Iss_avg8, img_gt); 
    [~,avg8_SSIMs(ii,3)] = comp_psnr_2(Ims_avg8, img_gt);
    [~,avg8_SSIMs(ii,4)] = comp_psnr_2(Ifd_avg8, img_gt);
    
    [~,avg16_SSIMs(ii,1)] = comp_psnr_2(img_avg16, img_gt);
    [~,avg16_SSIMs(ii,2)] = comp_psnr_2(Iss_avg16, img_gt); 
    [~,avg16_SSIMs(ii,3)] = comp_psnr_2(Ims_avg16, img_gt); 
    [~,avg16_SSIMs(ii,4)] = comp_psnr_2(Ifd_avg16, img_gt); 

%% Save images

        if ~isempty(save_dir)
            
            disp('Save Images ...')
            
            if ~exist(save_dir1,'dir')
                mkdir(save_dir1)
                mkdir(save_dir2)
                mkdir(save_dir4)
                mkdir(save_dir8)
                mkdir(save_dir16)
            end
            

            % Save output image
            ext_ = '.tif';
            
            name_str = [img_name postfix_Iss ext_];
            imwrite(Iss_avg1,fullfile(save_dir1,name_str))
            name_str = [img_name postfix_Ims ext_];
            imwrite(Ims_avg1,fullfile(save_dir1,name_str))
            name_str = [img_name postfix_Ifd ext_];
            imwrite(Ifd_avg1,fullfile(save_dir1,name_str))
            
            name_str = [img_name postfix_Iss ext_];
            imwrite(Iss_avg2,fullfile(save_dir2,name_str))
            name_str = [img_name postfix_Ims ext_];
            imwrite(Ims_avg2,fullfile(save_dir2,name_str))
            name_str = [img_name postfix_Ifd ext_];
            imwrite(Ifd_avg2,fullfile(save_dir2,name_str))
            
            name_str = [img_name postfix_Iss ext_];
            imwrite(Iss_avg4,fullfile(save_dir4,name_str))
            name_str = [img_name postfix_Ims ext_];
            imwrite(Ims_avg4,fullfile(save_dir4,name_str))
            name_str = [img_name postfix_Ifd ext_];
            imwrite(Ifd_avg4,fullfile(save_dir4,name_str))
            
            name_str = [img_name postfix_Iss ext_];
            imwrite(Iss_avg8,fullfile(save_dir8,name_str))
            name_str = [img_name postfix_Ims ext_];
            imwrite(Ims_avg8,fullfile(save_dir8,name_str))
            name_str = [img_name postfix_Ifd ext_];
            imwrite(Ifd_avg8,fullfile(save_dir8,name_str))
            
            name_str = [img_name postfix_Iss ext_];
            imwrite(Iss_avg16,fullfile(save_dir16,name_str))
            name_str = [img_name postfix_Ims ext_];
            imwrite(Ims_avg16,fullfile(save_dir16,name_str))
            name_str = [img_name postfix_Ifd ext_];
            imwrite(Ifd_avg16,fullfile(save_dir16,name_str))
            
            
            %{ 
            % Save ground-truth image
            
            to_save_fpath = fullfile(save_dir,[img_name ext_]);
            if ~exist(to_save_fpath,'file')
                imwrite(im/maxI,to_save_fpath)
            end
            %}
        end

end



%% 


% diary off



%%
% SAVE "avg1_TIMEs"
save(mat_fpath,...
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
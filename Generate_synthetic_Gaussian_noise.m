%{
Given a folder containing images, this script generates noisy version of
those images. For each image, a separate folder is created in the
destination directory (save_dir). Then, for each noise level, a subfolder
is created.

Example:
    Dataset folder contains: image1.png
    Input noise levels are indicated by "input_sigmas = [25];"
    The output directoy is as follows:
        output_directory\image1\sigma025\I1.mat
        output_directory\image1\sigma025\I7.mat
    where, "I1.mat" stores noisy image and "I7.mat" stores its
    ground-truth image.

%}



%%
clear
clc

data_path = './DATASETS/CSR_test_gray';

save_dir = './DATASETS/CSR_test_gray';



%% Destination directory

if ~exist(save_dir,'dir')
    mkdir(save_dir)
end



%% Noise levels

input_sigmas = [15,25,50];
% input_sigmas = 25;
N_sigmas = length(input_sigmas);
sigma_level = cell(N_sigmas,1);



%% list dataset directory to get file names

listing = dir(data_path);
N_images = length(listing) - 2;


%% Generate noisy images and store them in appropriate subfolders.

for i = 1:N_images
    image_name = listing(i+2).name;
    image_path = fullfile(data_path,image_name);
    
    [folder_path,image_name,ext] = fileparts(image_path);
    
    msg = '********* image %s (#%d/%d) *********\n\n';
    fprintf(msg,image_name,i,N_images); 
    
    
    for j = 1:N_sigmas
        
        sig = input_sigmas(j);
        
        fprintf('--------> noise %d (#%d/%d) \n\n',sig,j,N_sigmas);
        
        subfolder_name = sprintf('sigma%0.3d',sig);
        save_dir_noise = fullfile(save_dir,image_name,subfolder_name);
        
        if ~exist(save_dir_noise,'dir')
            mkdir(save_dir_noise)
        end 
        
        
        % Ground-truth image
        I7 = double(imread(image_path));
        
        
        % Generate NOISY image
        NOISE = randn(size(I7)) * sig;
        I1 = I7 + NOISE;
        
        save(fullfile(save_dir_noise,'I7.mat'),'I7') % GT
        save(fullfile(save_dir_noise,'I1.mat'),'I1') % NOISY
         
    end
    
end
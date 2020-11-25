%{

NOTE: REMOVE GSRC_NLP_Demo from path.


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
    'Ims','runtime.Ims','_%.3d_1_Ims';
    };
args.data_path = './RESULTS/blsgsm_dt12';
args.save_dir = './RESULTS/blsgsm_dt12_imgs'; % or ''

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
% close all
clc

addpath('./METHODS/myFunctions')
addpath('./METHODS/BLS-GSM')
addpath('./METHODS/BLS-GSM/Added_PyrTools')
addpath('./METHODS/BLS-GSM/denoising_subprograms')



data_folder = 'dt12'


data_path = fullfile('./DATASETS/',data_folder); %%%%%%%%%%%%%

save_dir = fullfile('./RESULTS/',['blsgsm_' data_folder]); %%%%%%%%%%%%% 


%% Generate folder name for each noise level

input_sigmas = [15,25,50]; %%%%%%%%%%%%%

N_sigmas = length(input_sigmas);




%% Reading all noisy and clean images in the dataset

listing = dir(data_path);
N_folders = length(listing) - 2;

for i = 1 % %%%%%%%%%%%%%

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
        
        [Ny,Nx] = size(imn);
        
        PS = ones(size(imn));	% power spectral density (in this case, flat, i.e., white noise)
        seed = 0;               % random seed

        % Pyramidal representation parameters
        Nsc = ceil(log2(min(Ny,Nx)) - 4);  % Number of scales (adapted to the image size)
        Nor = 3;				            % Number of orientations (for X-Y separable wavelets it can only be 3)
        repres1 = 'uw';                     % Type of pyramid (shift-invariant version of an orthogonal wavelet, in this case)
        repres2 = 'daub1';                  % Type of wavelet (daubechies wavelet, order 2N, for 'daubN'; in this case, 'Haar')
        
        % Model parameters (optimized: do not change them unless you are an advanced user with a deep understanding of the theory)
        blSize = [3 3];	    % n x n coefficient neighborhood of spatial neighbors within the same subband
                            % (n must be odd): 
        parent = 0;			% including or not (1/0) in the neighborhood a coefficient from the same spatial location
                            % and orientation as the central coefficient of the n x n neighborhood, but
                            % next coarser scale. Many times helps, but not always.
        boundary = 1;		% Boundary mirror extension, to avoid boundary artifacts 
        covariance = 1;     % Full covariance matrix (1) or only diagonal elements (0).
        optim = 1;          % Bayes Least Squares solution (1), or MAP-Wiener solution in two steps (0)

        % -------------------------
        % Uncomment the following 4 code lines for reproducing the results of our IEEE Trans. on Im. Proc., Nov. 2003 paper
        % This configuration is slower than the previous one, but it gives slightly better results (SNR)
        % on average for the test images "lena", "barbara", and "boats" used in the cited article.

        Nor = 8;                           % 8 orientations
        repres1 = 'fs';                    % Full Steerable Pyramid, 5 scales for 512x512
        repres2 = '';                      % Dummy parameter when using repres1 = 'fs'   
        parent = 1;                        % Include a parent in the neighborhood
        % -------------------------
        
        
        % Call the denoising function
        tic
        Ims = denoi_BLS_GSM(imn, sig, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed);
        runtime.Ims = toc;
        
%         imshow_eval_2(im,Ims,'BLS-GSM')


        [MM,NN] = size(im);
        ps_Ims = comp_psnr_2(im,Ims);
            
        fprintf('%9s\n','Ims')
        fprintf('%9.2f\n',ps_Ims)
        
        
%% Save results           
        
        %=================
        % Save MAT files %
        % =================
        
        save(outputs_path,...
            'im_path','imn_path',...            % clean and noisy image paths
            'sig',...                           % Gaussian noise level
            'Ims',...
            'runtime');  
        
    end
    
end

disp('FINISHED')
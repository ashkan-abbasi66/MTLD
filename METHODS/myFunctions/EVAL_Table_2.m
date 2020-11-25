%{ 
=== NOTE ===
Before using this script, 
 use "RENAME_FOLDERS" for "DATA_PATH" if needed

This function assumes all outputs are saved into a .mat file with name 
"outputs.mat".

Edit these lines: Search for "%%%%%%%%%%%%%".
============

*** There are three important parts in this file. 
1- Creating a table
2- Filling the table
2-1- inserting information

To create the table, three cell arrays are used: img_names, and cols_psnr [or
cols_ssim], and cols_time.


USAGE EXAMPLE:

clear 
clc
args.data_path = './RESULTS/TL_paper_extended_matFiles';
args.save_dir = './RESULTS/TL_paper_extended_images'; % or ''
args.image_variable_name = 'Iss'; %%%%%%%%%%%%%
args.runtime_variable_name = 'runtime.ksvd.Iss'; %%%%%%%%%%%%%
args.postfix = '_%.3d_1_ksvd'; % %%%%%%%%%%%%% _sigma_(alg. number)_(alg. name) 

[table_psnr,table_ssim,table_time] = EVAL_Table_2(args);

%}

function [table_psnr,table_ssim,table_time] = EVAL_Table_2(args)

% ==================================
% Get input arguments
% ==================================
data_path = args.data_path;
if isfield(args,'save_dir')
    save_dir = args.save_dir;
else
    save_dir = '';
end
image_var_name = args.image_variable_name;
if isfield(args,'runtime_variable_name')
    time_var_name = args.runtime_variable_name;
else
    time_var_name = '';
end
postfix = args.postfix;


fnames = dir(data_path);
N_folders = length(fnames) - 2;

% ==================================
% Create a Table
% ==================================


% Get Distinct Image Names %

img_names = {'Image Name'};

for i = 1:N_folders
    image_name = fnames(i+2).name;   
    img_names{i+1,1} = image_name;    
end


% Get Noise Levels%

noise_levels = dir(fullfile(data_path,img_names{2}));
N_sig = length(noise_levels) - 2;

cols_psnr = cell(N_folders + 1,N_sig);
for i = 1:N_sig
    str = noise_levels(i+2).name;
    str_ = str(strfind(str,'a')+1:end);
    sig = str2double(str_);
    cols_psnr{1,i} = sig;
end

cols_ssim = cols_psnr;

cols_time = cols_psnr;



% ==================================
% Fill the Table
% ==================================

for i = 1:N_folders
    
    img_name = img_names{i+1,1};
    
    fprintf('====> IMAGE %s (#%d/%d) \n',img_name,i,N_folders);
    
    for j = 1:N_sig
        
        sigma_str = noise_levels(j+2).name;
        str = sigma_str;
        str_ = str(strfind(str,'a')+1:end);
        sig = str2double(str_);
        
        fprintf('      sigma %d (#%d/%d) \n',sig,j,N_sig);
        
        
        
        % Load files
        fpath = fullfile(data_path,img_name,sigma_str,'outputs.mat');
        load(fpath);
        
        load(imn_path)
        load(im_path)
        
        imn = I1(:,:,1); %%%%%%%%%%%%%
        im = I7(:,:,1); %%%%%%%%%%%%%
        
        clear I1 I7
        
        
        
        % find column number
        for k = 1:size(cols_psnr,2)
            if sig == cols_psnr{1,k}
                sig_ind = k;
                break;
            end
        end
    
        %find row number
        for k = 1:length(img_names)
            if strcmpi(img_name,img_names{k,1})
                img_ind = k;
                break;
            end
        end
        
        
        
        % ==================================
        % insert information at location (img_ind,sig_ind) %
        % ==================================
        if ~strcmpi(image_var_name,'imn')
            result_img = eval(image_var_name); %%%%%%%%%%%%% Step (1/3)
        else
            result_img = imn;
        end

        [psnr_,ssim_] = comp_psnr_2(im,result_img(1:size(im,1),1:size(im,2)));
        cols_psnr{img_ind,sig_ind} = psnr_;
        cols_ssim{img_ind,sig_ind} = ssim_;
        if ~isempty(time_var_name)
            cols_time{img_ind,sig_ind} = eval(time_var_name); %%%%%%%%%%%%% Step (2/3)
        end
        
        if ~isempty(save_dir)
%             postfix = '_%.3d_6_mx_swt_121'; %%%%%%%%%%%%% Step (3/3)
            
            if ~exist(save_dir,'dir')
                mkdir(save_dir)
            end

            % Save output image
            maxI = 255;
            ext_ = '.tif';
            name_str = [img_name postfix ext_];
            to_save_fpath = fullfile(save_dir,sprintf(name_str,sig));
            imwrite(result_img/maxI,to_save_fpath)
            %imwrite(uint8(result_img),to_save_fpath,'tif');

            % Save ground-truth image
            to_save_fpath = fullfile(save_dir,[img_name ext_]);
            if ~exist(to_save_fpath,'file')
                imwrite(im/maxI,to_save_fpath)
                %imwrite(uint8(im),to_save_fpath,'tif');
            end
        end
        
    end
end

table_psnr = table(img_names,cols_psnr);
table_ssim = table(img_names,cols_ssim);
table_time = table(img_names,cols_time);

end
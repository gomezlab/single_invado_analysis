function make_comparison_frames(cfg_file)

[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

addpath(genpath(path_folders));

if (exist('comp_base_out_folder','var'))
    if (not(exist(comp_base_out_folder,'dir')))
        mkdir(comp_base_out_folder);
    end
end

extr_val_fullfile = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_fullfile,'file'))
    min_max = load(extr_val_fullfile);
end

i_count = size(imfinfo(original_i_file),2);
pad_length = num2str(length(num2str(i_count)));
i_size = size(imread(original_i_file,1));
pix_count = i_size(1)*i_size(2);

for i = 1:i_count
    padded_i_num = sprintf(['%0',pad_length,'d'],i);
    this_data_folder = fullfile(comp_base_in_folder,padded_i_num);
    if (not(exist(this_data_folder,'dir')))
        continue;
    end
    if (i > 1), continue, end

    comp_image_nums = load(fullfile(comp_base_in_folder,padded_i_num,comp_image_nums_file))
    comp_starting_ad = load(fullfile(comp_base_in_folder,padded_i_num,comp_starting_adh_file));
    comp_ending_ad = load(fullfile(comp_base_in_folder,padded_i_num,comp_ending_adh_file));
    
    padded_i_num_1 = sprintf(['%0',pad_length,'d'],comp_image_nums(1));
    padded_i_num_2 = sprintf(['%0',pad_length,'d'],comp_image_nums(2));
    
    i1.orig = normalize_grayscale_image(imread(original_i_file,comp_image_nums(1)),min_max(1),min_max(2));
    i2.orig = normalize_grayscale_image(imread(original_i_file,comp_image_nums(2)),min_max(1),min_max(2));
    
    i1.cell_edge = bwperim(imread(fullfile(I_folder_2,padded_i_num_1,edge_filename)));
    i2.cell_edge = bwperim(imread(fullfile(I_folder_2,padded_i_num_2,edge_filename)));

    i1.adhesions = logical(imread(fullfile(I_folder_2,padded_i_num_1,t_filtered_file)));
    i1.label_ad  = bwlabel(i1.adhesions);
    i2.adhesions = logical(imread(fullfile(I_folder_2,padded_i_num_2,t_filtered_file)));
    i2.label_ad  = bwlabel(i2.adhesions);
    
    i1.adh_of_interest = zeros(size(i1.adhesions,1),size(i1.adhesions,2));
    i2.adh_of_interest = zeros(size(i2.adhesions,1),size(i2.adhesions,2));
    for j = 1:size(comp_starting_ad,2)
        i1.adh_of_interest(i1.label_ad == comp_starting_ad(j)) = 1;
        i2.adh_of_interest(i2.label_ad == comp_ending_ad(j))   = 1;
    end
    
     
%     
%     for j = 1:2
%         for k = 1:size(comp_image_nums(i))
% 
%             i_num = comp_image_nums(i,j);
% 
%             orig_i = normalize_grayscale_image(imread(original_i_file,i),min_max(1),min_max(2));
% 
%             adh_i = imread(fullfile(I_folder_1,i_num,t_filtered));
% 
%             I_lab = bwlabel(adh_i);
% 
%             highlighted_all = cat(3,orig_i,orig_i,orig_i);
%         end
%     end
end
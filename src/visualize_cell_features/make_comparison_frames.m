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

for i = 1:i_count
    padded_i_num = sprintf(['%0',pad_length,'d'],i);
    image_prob_folder = fullfile(comp_base_in_folder,padded_i_num);
    if (not(exist(image_prob_folder,'dir')))
        continue;
    end
    
    disp(padded_i_num)
    
    for j = 0:4000
        this_data_folder = fullfile(image_prob_folder,num2str(j));
        if (not(exist(this_data_folder,'dir'))), break, end
        
        comp_image_nums = load(fullfile(this_data_folder,comp_image_nums_file));
        comp_starting_ad = load(fullfile(this_data_folder,comp_starting_adh_file)) + 1;
        comp_ending_ad = load(fullfile(this_data_folder,comp_ending_adh_file)) + 1;
        comp_winning_ad = load(fullfile(this_data_folder,comp_winning_adh_file)) + 1;

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

        i1.adh_of_interest = zeros(i_size);
        i2.adh_of_interest = zeros(i_size);
        i1.highlighted = create_highlighted_image(i1.orig,i1.adhesions);
        i2.highlighted = create_highlighted_image(i2.orig,i2.adhesions);

        color_map = jet(size(comp_starting_ad,2)+1);
        winning_color = zeros(3,1);
        for k = 1:size(comp_starting_ad,2)
            if (comp_starting_ad(k) == comp_winning_ad)
                winning_color = color_map(k,:);
            end
        end
        
        for k = 1:size(comp_starting_ad,2)
            temp = zeros(i_size);
            temp(i1.label_ad == comp_starting_ad(k)) = 1;
            i1.adh_of_interest(i1.label_ad == comp_starting_ad(k)) = 1;
            i1.highlighted = create_highlighted_image(i1.highlighted,temp,'color',color_map(k,:));
            
            temp = zeros(i_size);
            temp(i2.label_ad == comp_ending_ad(k)) = 1;
            i2.adh_of_interest(i2.label_ad == comp_ending_ad(k))   = 1;
            i2.highlighted = create_highlighted_image(i2.highlighted,temp,'color',winning_color);
        end
        
        [min_x,min_y,max_x,max_y] = find_binary_bounding_box(i1.adh_of_interest + i2.adh_of_interest);
        
        min_x = min_x - image_padding_min;
        min_y = min_y - image_padding_min;
        max_x = max_x + image_padding_min;
        max_y = max_y + image_padding_min;
        
        if (min_x <= 0), min_x = 1; end
        if (min_y <= 0), min_x = 1; end
        if (max_x > i_size(2)), max_x = i_size(2); end
        if (max_y > i_size(1)), max_y = i_size(1); end
        
        if (not(exist(fullfile(comp_base_out_folder,padded_i_num),'dir'))) 
            mkdir(fullfile(comp_base_out_folder,padded_i_num));
        end
        
        i1.highlight_trim = i1.highlighted(min_y:max_y,min_x:max_x,1:3);
        i2.highlight_trim = i2.highlighted(min_y:max_y,min_x:max_x,1:3);
        trim_size = size(i1.highlight_trim);
        trim_size(2) = round(trim_size(2)*0.05);
        
        frame = [i1.highlight_trim,0.5*ones(trim_size),i2.highlight_trim];
        imwrite(frame,fullfile(comp_base_out_folder,padded_i_num,[num2str(j),'.png']));        
    end
end
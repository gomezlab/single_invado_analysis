function make_single_ad_frames(cfg_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'MAKE_SINGLE_AD_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

addpath(genpath(path_folders));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect General Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_image_num = find_max_image_num(I_folder);
folder_char_length = length(num2str(max_image_num));
i_size = size(imread(fullfile(I_folder,num2str(max_image_num),focal_image)));

tracking_seq = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Bounding Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bounding_matrix = [Inf*ones(size(tracking_seq,1),1), Inf*ones(size(tracking_seq,1),1), -Inf*ones(size(tracking_seq,1),1), -Inf*ones(size(tracking_seq,1),1)];

i_seen = 0;

for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);
    
    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;

    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));

    bounds = regionprops(ad_label,'BoundingBox');
    
    for i = 1:size(tracking_seq,1)
        tracking_row = tracking_seq(i,:);
        if (tracking_row(i_seen) <= 0), continue; end
    
        ad_num = tracking_row(i_seen);
        
        corners = [bounds(ad_num).BoundingBox(1), bounds(ad_num).BoundingBox(2)];
        corners = [corners, corners + bounds(ad_num).BoundingBox(3:4)];
        
        if (corners(1) < bounding_matrix(i,1)), bounding_matrix(i,1) = corners(1); end
        if (corners(2) < bounding_matrix(i,2)), bounding_matrix(i,2) = corners(2); end
        if (corners(3) > bounding_matrix(i,3)), bounding_matrix(i,3) = corners(3); end
        if (corners(4) > bounding_matrix(i,4)), bounding_matrix(i,4) = corners(4); end
    end
end

bounding_matrix(1,1) = -1;
bounding_matrix(:,1:2) = bounding_matrix(:,1:2) - image_padding_min;
bounding_matrix(bounding_matrix(:,1) <= 0,1) = 1;
bounding_matrix(find(bounding_matrix(:,2) <= 0),2) = 1;

bounding_matrix(:,3:4) = bounding_matrix(:,3:4) + image_padding_min;
bounding_matrix(find(bounding_matrix(:,3) > i_size(2)),3) = i_size(2);
bounding_matrix(find(bounding_matrix(:,4) > i_size(1)),4) = i_size(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Single Ad Image Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_seen = 0;
i_visualized_mat = zeros(1,size(tracking_seq,1));

for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;
    
    %Gather and scale the input adhesion image
    orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
    scale_factor = double(intmax(class(orig_i)));
    orig_i = double(orig_i)/scale_factor;

    %Gather and process the ad label image    
    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
    ad_label_perim = zeros(i_size);
    for i = 1:max(ad_label(:))
        assert(any(any(ad_label == i)), 'Error: can''t find ad number %d, in adhesion image number %d',i,padded_i_num);
        
        this_ad = zeros(i_size);
        this_ad(ad_label == i) = 1;
        ad_label_perim(bwperim(this_ad)) = i;
    end
    
    %Gather the cell edge image if available
    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
    end
    
    for i = 1:size(tracking_seq,1)
        tracking_row = tracking_seq(i,:);
        if (tracking_row(i_seen) <= 0), continue; end
        
        i_visualized_mat(i) = i_visualized_mat(i) + 1;
        
        ad_num = tracking_row(i_seen);
        this_ad = zeros(size(ad_label));
        this_ad(ad_label == ad_num) = 1;

        not_this_ad = xor(im2bw(ad_label,0),this_ad);
        assert(sum(sum(not_this_ad)) + sum(sum(this_ad)) == sum(sum(im2bw(ad_label,0))))
    end
    if (mod(j,10) == 0), disp(j); end
end


1;

%     padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);
% 
%     if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end
% 
%     i_seen = i_seen + 1;
%     padded_i_seen = sprintf(['%0',num2str(folder_char_length),'d'],i_seen);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Gather and scale the input adhesion image
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
%     scale_factor = double(intmax(class(orig_i)));
%     orig_i = double(orig_i)/scale_factor;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Gather the adhesion label image and perimeters
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
%     ad_nums = tracking_seq(tracking_seq(:,i_seen) > 0,i_seen);
%     for j = 1:length(ad_nums)
%         this_num = ad_nums(j);
%         assert(any(any(ad_label == this_num)), 'Error: can''t find ad num %d in image number %d.',this_num,padded_i_num)
%     end
% 
%     ad_label_perim = zeros(size(orig_i,1),size(orig_i,2));
%     for j = 1:length(ad_nums)
%         this_num = ad_nums(j);
%         this_ad = zeros(size(orig_i,1),size(orig_i,2));
%         this_ad(ad_label == this_num) = 1;
%         ad_label_perim(bwperim(this_ad)) = this_num;
%     end
% 
%     if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
%         cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Other images
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %Build the unique lineage highlighted image
%     cmap_nums = lineage_to_cmap(tracking_seq(:,i_seen) > 0);
%     assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',padded_i_num);
%     this_cmap = zeros(max(ad_label_perim(:)),3);
%     this_cmap(ad_nums,:) = lineage_cmap(cmap_nums,:);
% 
%     highlighted_all = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);
% 
%     %Build the birth time highlighted image
%     cmap_nums = birth_time_to_cmap(tracking_seq(:,i_seen) > 0);
%     assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in birth time image %d',padded_i_num);
%     this_cmap = zeros(max(ad_label_perim(:)),3);
%     this_cmap(ad_nums,:) = time_cmap(cmap_nums,:);
% 
%     highlighted_time = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);
% 
%     if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
%         cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
%         edge_image_ad = create_highlighted_image(edge_image_ad,cell_edge,'color_map',edge_cmap(i_seen,:));
%     end
%     edge_image_ad = create_highlighted_image(edge_image_ad,im2bw(ad_label_perim,0),'color_map',edge_cmap(i_seen,:));
% 
%     orig_i = orig_i(b_box(2):b_box(4), b_box(1):b_box(3));
%     highlighted_all = highlighted_all(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
%     highlighted_time = highlighted_time(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
%     edge_image_ad_bounded = edge_image_ad(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
% 
%     spacer = 0.5*ones(size(orig_i,1),1,3);
% 
%     if(i_p.Results.debug), disp(i_seen); end

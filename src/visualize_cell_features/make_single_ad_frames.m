function make_single_ad_frames(cfg_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'MAKE_SINGLE_AD_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});
if (i_p.Results.debug == 1), profile off; profile on; end

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
        corners = [corners, corners + bounds(ad_num).BoundingBox(3:4)]; %#ok<AGROW>

        if (corners(1) < bounding_matrix(i,1)), bounding_matrix(i,1) = corners(1); end
        if (corners(2) < bounding_matrix(i,2)), bounding_matrix(i,2) = corners(2); end
        if (corners(3) > bounding_matrix(i,3)), bounding_matrix(i,3) = corners(3); end
        if (corners(4) > bounding_matrix(i,4)), bounding_matrix(i,4) = corners(4); end
    end
end

bounding_matrix(:,1:2) = bounding_matrix(:,1:2) - image_padding_min;
bounding_matrix(bounding_matrix(:,1) <= 0,1) = 1;
bounding_matrix(bounding_matrix(:,2) <= 0,2) = 1;
bounding_matrix(:,1:2) = floor(bounding_matrix(:,1:2));

bounding_matrix(:,3:4) = bounding_matrix(:,3:4) + image_padding_min;
bounding_matrix(bounding_matrix(:,3) > i_size(2),3) = i_size(2);
bounding_matrix(bounding_matrix(:,4) > i_size(1),4) = i_size(1);
bounding_matrix(:,3:4) = ceil(bounding_matrix(:,3:4));

lineage_lengths = zeros(1,size(tracking_seq,1));
for i = 1:size(tracking_seq,1)
    lineage_lengths(i) = sum(tracking_seq(i,:) > 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Single Ad Image Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_seen = 0;

all_images = cell(size(tracking_seq,1), max_image_num);
for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;
    surrounding_cols = zeros(size(tracking_seq,1),3);
    surrounding_cols(:,2) = tracking_seq(:,i_seen);
    try, surrounding_cols(:,1) = tracking_seq(:,i_seen - 1); end
    try, surrounding_cols(:,3) = tracking_seq(:,i_seen + 1); end
    
    if (all(all(surrounding_cols <= 0))), continue; end
    
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
        surrounding_entries = [0, tracking_row(i_seen), 0];
        try, surrounding_entries(1) = tracking_row(i_seen - 1); end
        try, surrounding_entries(3) = tracking_row(i_seen + 1); end
        
        if (all(surrounding_entries <= 0)), continue; end

        ad_num = tracking_row(i_seen);
        if (ad_num <= 0); ad_num = -Inf; end
        bounded_ad_label_perim = ad_label_perim(bounding_matrix(i,2):bounding_matrix(i,4), bounding_matrix(i,1):bounding_matrix(i,3));

        this_ad = zeros(size(bounded_ad_label_perim));
        this_ad(bounded_ad_label_perim == ad_num) = 1;

        not_this_ad = xor(im2bw(bounded_ad_label_perim,0),this_ad);
        assert(sum(sum(not_this_ad)) + sum(sum(this_ad)) == sum(sum(im2bw(bounded_ad_label_perim,0))))

        highlighted_image = orig_i(bounding_matrix(i,2):bounding_matrix(i,4), bounding_matrix(i,1):bounding_matrix(i,3));
        highlighted_image = create_highlighted_image(highlighted_image,this_ad,'color_map',[0,1,0],'mix_percent',0.5);
        highlighted_image = create_highlighted_image(highlighted_image,not_this_ad,'color_map',[0,0,1],'mix_percent',0.5);
        if (exist('cell_edge','var'))
            bounded_edge = cell_edge(bounding_matrix(i,2):bounding_matrix(i,4), bounding_matrix(i,1):bounding_matrix(i,3));
            highlighted_image = create_highlighted_image(highlighted_image,bounded_edge,'color_map',[1,0,0],'mix_percent',0.5);
        end

        all_images{i,j} = highlighted_image;
    end
    if (mod(j,10) == 0 && i_p.Results.debug), disp(j); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Ad Sequences in Montage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(all_images,1)
    this_set = all_images(i,1:end);
    notempty = zeros(size(this_set));
    for j = 1:size(this_set,2)
        notempty(j) = not(isempty(this_set{j}));
    end
    if (sum(notempty) == 0), continue; end

    total_images = find(notempty,1,'last') - find(notempty,1,'first') + 1;

    images_per_side = ones(1,2);
    for j = 1:(ceil(sqrt(total_images)) - 1)
        if (total_images <= images_per_side(1)*images_per_side(2)), continue; end

        if (total_images <= j * (j + 1))
            images_per_side = [j, j+1];
            continue;
        else
            images_per_side = [j + 1, j + 1];
        end
    end
    assert(images_per_side(1)*images_per_side(2) >= total_images, 'Error; images per side not large enough');
    image_size = size(this_set{1,find(notempty,1,'first')});

    montage = 0.5*ones(image_size(1)*images_per_side(1)+images_per_side(1)-1, image_size(2)*images_per_side(2)+images_per_side(2) - 1, 3);
    for j = 1:images_per_side(1)
        for k = 1:images_per_side(2)
            i_index = find(notempty,1,'first') + (j-1)*images_per_side(2) + (k-1);

            if (i_index > length(notempty)), continue; end
            if (not(notempty(i_index))), continue; end

            montage((j-1)*image_size(1) + j:(j)*image_size(1) + j - 1, (k-1)*image_size(2) + k:(k)*image_size(2) + k - 1, 1:3) = this_set{1,i_index};
        end
    end
    if (not(exist(fullfile(out_path,'single_ad'),'dir'))), mkdir(fullfile(out_path,'single_ad')); end

    padded_num = sprintf(['%0',num2str(length(num2str(size(all_images,1)))),'d'],i);
    imwrite(montage,fullfile(out_path,'single_ad', ['montage_',padded_num, '.png']));
end

profile off;
if (i_p.Results.debug), profile viewer; end
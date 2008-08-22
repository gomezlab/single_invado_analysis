function make_movie_frames(cfg_file,varargin)
%MAKE_MOVIE_FRAMES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'MAKE_MOVIE_FRAMES';

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
I_folder_dir = dir(I_folder);
max_image_num = -Inf;
folder_char_length = 0;
i_size = 0;
for i = 1:size(I_folder_dir)
    if (strcmp(I_folder_dir(i).name,'.') || strcmp(I_folder_dir(i).name,'..') || not(I_folder_dir(i).isdir)), continue; end

    if (str2num(I_folder_dir(i).name) > max_image_num), max_image_num = str2num(I_folder_dir(i).name); end

    folder_char_length = length(I_folder_dir(i).name);

    poss_sample_image =fullfile(I_folder,sprintf(['%0',num2str(folder_char_length),'d'],i),focal_image);
    if (exist(poss_sample_image,'file') && length(i_size) == 1)
        i_size = size(imread(poss_sample_image));
    end
end

tracking_seq = load(tracking_seq_file) + 1;

ghost_frames_count = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist(bounding_box_file,'file'))
    b_box = load(bounding_box_file);
else
    b_box = find_time_series_bbox(I_folder);
    csvwrite(bounding_box_file,b_box);
end

b_box(1:2) = b_box(1:2) - image_padding_min;
b_box(3:4) = b_box(3:4) + image_padding_min;
if (b_box(1) <= 0), b_box(1) = 1; end
if (b_box(2) <= 0), b_box(2) = 1; end
if (b_box(3) > i_size(2)), b_box(3) = i_size(2); end
if (b_box(4) > i_size(1)), b_box(4) = i_size(1); end

edge_cmap = jet(size(tracking_seq,2));
%define the edge image here because the edge image will be added to each
%image loop, so the image should be global
edge_image_ad = ones(i_size(1),i_size(2),3);

max_live_adhesions = find_max_live_adhesions(tracking_seq);

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);

time_cmap = jet(size(tracking_seq,2));
birth_time_to_cmap = zeros(size(tracking_seq,1),1);

i_seen = 0;

for i = 1:max_image_num
    if (i_seen + 1 > size(tracking_seq,2))
        continue;
    end

    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;
    padded_i_seen = sprintf(['%0',num2str(folder_char_length),'d'],i_seen);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather and scale the input adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
    scale_factor = double(intmax(class(orig_i)));
    orig_i = double(orig_i)/scale_factor;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
    ad_nums = tracking_seq(tracking_seq(:,i_seen) > 0,i_seen);
    for j = 1:length(ad_nums)
        this_num = ad_nums(j);
        assert(any(any(ad_label == this_num)), 'Error: can''t find ad num %d in image number %d.',this_num,padded_i_num)
    end

    ad_label_perim = zeros(size(orig_i,1),size(orig_i,2));
    for j = 1:length(ad_nums)
        this_num = ad_nums(j);
        this_ad = zeros(size(orig_i,1),size(orig_i,2));
        this_ad(ad_label == this_num) = 1;
        ad_label_perim(bwperim(this_ad)) = this_num;
    end

    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the matrices translating number to colormap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(tracking_seq,1)
        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_seen) <= 0), continue; end

        %Add the current adhesion to the time dependent color map if no
        %number is currently defined
        if (birth_time_to_cmap(j) == 0), birth_time_to_cmap(j) = i_seen; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_seq(:,i_seen) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);

            taken_nums = zeros(1,max_live_adhesions);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);

            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
            catch
                assert(isempty(find(taken_dists == max(taken_dists),1,'first')), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end

    %Make sure all the live adhesions have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning unique color codes');
    assert(all(birth_time_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning birth time color codes');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Adhesion Ghost Images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Save the label matrices
    if (exist('labels','var'))
        frame_size_count = size(labels,2);
        if (frame_size_count >= ghost_frames_count), frame_size_count = ghost_frames_count - 1; end
        for j = frame_size_count:-1:1
            labels(j+1).ad_perim = labels(j).ad_perim;
            labels(j+1).ad_filled = labels(j).ad_filled;
        end
        labels(1).ad_perim = ad_label_perim;
        labels(1).ad_filled = ad_label;
    else
        labels(1).ad_perim = ad_label_perim;
        labels(1).ad_filled = ad_label;
    end

    %Draw the ghost images
    if (i_seen == size(tracking_seq,2))
        highlighted_ghost_unique = zeros(size(orig_i));
        highlighted_ghost_time = zeros(size(orig_i));

        highlighted_ghost_unique_filled = zeros(size(orig_i));
        highlighted_ghost_time_filled = zeros(size(orig_i));
        for m=size(labels,2):-1:1
            this_i_num = i_seen - m + 1;
            this_ad_perim = labels(m).ad_perim;
            this_ad_filled = labels(m).ad_filled;

            these_ad_nums = tracking_seq(tracking_seq(:,this_i_num) > 0,this_i_num);

            mix_percent = (size(labels,2) - m + 1)/size(labels,2);

            %Unique colored adhesion lineage image drawing
            cmap_nums = lineage_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in unique color ghost image creation %d',this_i_num);
            this_cmap = zeros(max(this_ad_perim(:)),3);
            for j=1:length(cmap_nums)
                this_cmap(these_ad_nums(j),:) = lineage_cmap(cmap_nums(j),:);
            end
            highlighted_ghost_unique = create_highlighted_image(highlighted_ghost_unique,this_ad_perim,'color_map',this_cmap,'mix_percent',mix_percent);
            highlighted_ghost_unique_filled = create_highlighted_image(highlighted_ghost_unique_filled,this_ad_filled,'color_map',this_cmap,'mix_percent',mix_percent);

            %Birth time colored adhesion lineage image drawing
            cmap_nums = birth_time_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in birth time color ghost image creation %d',this_i_num);
            this_cmap = zeros(max(this_ad_perim(:)),3);
            for j=1:length(cmap_nums)
                this_cmap(these_ad_nums(j),:) = time_cmap(cmap_nums(j),:);
            end
            highlighted_ghost_time = create_highlighted_image(highlighted_ghost_time,this_ad_perim,'color_map',this_cmap,'mix_percent',mix_percent);
            highlighted_ghost_time_filled = create_highlighted_image(highlighted_ghost_time_filled,this_ad_filled,'color_map',this_cmap,'mix_percent',mix_percent);
        end
        highlighted_ghost_unique = highlighted_ghost_unique(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_time = highlighted_ghost_time(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_unique_filled = highlighted_ghost_unique_filled(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_time_filled = highlighted_ghost_time_filled(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        if (exist('pixel_size','var'))
            highlighted_ghost_unique = draw_scale_bar(highlighted_ghost_unique,pixel_size);
            highlighted_ghost_time = draw_scale_bar(highlighted_ghost_time,pixel_size);
            highlighted_ghost_unique_filled = draw_scale_bar(highlighted_ghost_unique_filled,pixel_size);
            highlighted_ghost_time_filled = draw_scale_bar(highlighted_ghost_time_filled,pixel_size);
        end
        imwrite(highlighted_ghost_unique,fullfile(out_path,'ghost_uni.png'));
        imwrite(highlighted_ghost_time,fullfile(out_path,'ghost_time.png'));
        imwrite(highlighted_ghost_unique_filled,fullfile(out_path,'ghost_uni_filled.png'));
        imwrite(highlighted_ghost_time_filled,fullfile(out_path,'ghost_time_filled.png'));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Other images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i_seen) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',padded_i_num);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = lineage_cmap(cmap_nums,:);

    highlighted_all = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);

    %Build the birth time highlighted image
    cmap_nums = birth_time_to_cmap(tracking_seq(:,i_seen) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in birth time image %d',padded_i_num);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = time_cmap(cmap_nums,:);

    highlighted_time = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);

    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
        edge_image_ad = create_highlighted_image(edge_image_ad,cell_edge,'color_map',edge_cmap(i_seen,:));
    end
    edge_image_ad = create_highlighted_image(edge_image_ad,im2bw(ad_label_perim,0),'color_map',edge_cmap(i_seen,:));

    orig_i = orig_i(b_box(2):b_box(4), b_box(1):b_box(3));
    highlighted_all = highlighted_all(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
    highlighted_time = highlighted_time(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
    edge_image_ad_bounded = edge_image_ad(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);

    spacer = 0.5*ones(size(orig_i,1),1,3);

    frame = cell(1,3);
    frame{1} = [edge_image_ad_bounded,spacer,highlighted_all];
    frame{2} = [cat(3,orig_i,orig_i,orig_i),spacer,highlighted_all];
    frame{3} = [edge_image_ad_bounded,spacer,highlighted_time];

    if (exist('pixel_size','var'))
        for j = 1:size(frame,2)
            frame{j} = draw_scale_bar(frame{j},pixel_size);
        end
    end

    if (exist('out_path','var'))
        for j = 1:length(out_prefix)
            output_filename = fullfile(out_path,out_prefix{1,j},[padded_i_seen,'.png']);
            fullpath = fileparts(output_filename);
            if (not(exist(fullpath,'dir')))
                mkdir(fullpath);
            end
            imwrite(frame{j},output_filename);
        end
    end

    if(i_p.Results.debug), disp(i_seen); end
end

profile off;
if (i_p.Results.debug), profile viewer; end

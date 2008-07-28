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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist(bounding_box_file,'file'))
    bounding_box = load(bounding_box_file);
else
    bounding_box = find_time_series_bbox(I_folder);
    csvwrite(bounding_box_file,bounding_box);
end

bounding_box(1:2) = bounding_box(1:2) - image_padding_min;
bounding_box(3:4) = bounding_box(3:4) + image_padding_min;
if (bounding_box(1) <= 0), bounding_box(1) = 1; end
if (bounding_box(2) <= 0), bounding_box(2) = 1; end
if (bounding_box(3) > i_size(2)), bounding_box(3) = i_size(2); end
if (bounding_box(4) > i_size(1)), bounding_box(4) = i_size(1); end

edge_cmap = jet(size(tracking_seq,2));
%define the edge image here because the old status of the edge image should
%be preserved
edge_image_ad = ones(i_size(1),i_size(2),3);

max_live_adhesions = find_max_live_adhesions(tracking_seq);

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);

time_cmap = jet(size(tracking_seq,2));
birth_time_to_cmap = zeros(size(tracking_seq,1),1);

i_seen = 0;

old_frames_count = 1000;

for i = 1:max_image_num
    if (i_seen + 1 > size(tracking_seq,2))
        continue;
    end

    i_seen = i_seen + 1;

    if (i_p.Results.debug && i_seen > 1000); continue; end

    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);
    padded_i_seen = sprintf(['%0',num2str(folder_char_length),'d'],i_seen);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end
    
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
    
    %some tracking matrices do not include all the adhesions, but much
    %of the subsequent code assumes that all ad numbers between 1 and the
    %number of adhesions are present, so relabel the adhesions where
    %numbers are missing in the tracking matrix column
    if (length(ad_nums) ~= length(unique(ad_label)))
        ad_label_temp = zeros(size(ad_label));
        for j=1:length(ad_nums)
            ad_label_temp(ad_label == ad_nums(j)) = j;
        end
        ad_label = ad_label_temp;
        ad_nums = unique(ad_label);
        if (ad_nums(1) == 0), ad_nums = ad_nums(2:end); end
    end
    
    ad_label_perim = zeros(size(orig_i,1),size(orig_i,2));
    for j=1:max(ad_label(:))
        this_ad = zeros(size(orig_i,1),size(orig_i,2));
        this_ad(ad_label == j) = 1;
        ad_label_perim(bwperim(this_ad)) = j;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the color map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    search_ascending = 1;
    for j = 1:size(tracking_seq,1)

        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_seen) <= 0), continue; end

        %Add the current adhesion to the time dependent color map if no
        %number is currently defined
        if (birth_time_to_cmap(j) == 0), birth_time_to_cmap(j) = i_seen; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            unique_used_colors = unique(lineage_to_cmap(tracking_seq(:,i_seen) > 0));

            if (search_ascending)
                search_seq = 1:max_live_adhesions;
                search_ascending = 0;
            else
                search_seq = max_live_adhesions:-1:1;
                search_ascending = 1;
            end

            poss_color = 0;
            for k = search_seq
                if (not(any(find(unique_used_colors == k))))
                    poss_color = k;
                    break;
                end
            end

            if (not(poss_color))
                warning('Matlab:adColor',['Could not find a unique color for adhesion lineage #: ', num2str(j)])
                lineage_to_cmap(j) = 1;
            else
                lineage_to_cmap(j) = poss_color;
            end
        end
    end

    %Make sure all the live adhesions have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning unique color codes');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Adhesion Ghost Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Save the label matrix
    if (exist('label_frames','var'))
        frame_size_count = size(label_frames,2);
        if (frame_size_count > old_frames_count), frame_size_count = old_frames_count; end
        for j = frame_size_count:-1:1
            label_frames{j + 1} = label_frames{j};
        end
        label_frames{1} = ad_label_perim;
    else
        label_frames{1} = ad_label_perim;
    end
    
    %Draw the ghost images
    if (i_seen == size(tracking_seq,2))
        highlighted_ghost_all = zeros(size(orig_i));
        highlighted_ghost_time = zeros(size(orig_i));
        for m=size(label_frames,2):-1:1
            
            this_i_num = i_seen - m + 1;
            
            these_ad_nums = tracking_seq(tracking_seq(:,this_i_num) > 0,this_i_num);
            
            labels = label_frames{m};

            mix_percent = (size(label_frames,2) - m + 1)/size(label_frames,2);

            cmap_nums = lineage_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in unique lineage numbers image %d',i);
            this_cmap = zeros(length(cmap_nums),3);
            for j=1:length(cmap_nums) 
                this_cmap(these_ad_nums(j),:) = lineage_cmap(cmap_nums(j),:);
            end
            highlighted_ghost_all = create_highlighted_image(highlighted_ghost_all,labels,'color_map',this_cmap,'mix_percent',mix_percent);
        
            cmap_nums = birth_time_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in unique lineage numbers image %d',i);
            this_cmap = zeros(length(cmap_nums),3);
            for j=1:length(cmap_nums) 
                this_cmap(these_ad_nums(j),:) = time_cmap(cmap_nums(j),:);
            end
            highlighted_ghost_time = create_highlighted_image(highlighted_ghost_time,labels,'color_map',this_cmap,'mix_percent',mix_percent);
        end
        imwrite(highlighted_ghost_all,fullfile(out_path,[num2str(i_seen),'_all','.png']));
        imwrite(highlighted_ghost_time,fullfile(out_path,[num2str(i_seen),'_time','.png']));        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Adhesion Ghost Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i_seen) > 0);

    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in unique lineage numbers image %d',i);

    all_cmap = zeros(length(cmap_nums),3);
    for j=1:length(cmap_nums)
        all_cmap(ad_nums(j),:) = lineage_cmap(cmap_nums(j),:);
    end
    highlighted_all = create_highlighted_image(orig_i,ad_label_perim,'color_map',all_cmap);
    
    %Build the birth time highlighted image
    cmap_nums = birth_time_to_cmap(tracking_seq(:,i_seen) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in birth time image %d',i);
    time_cmap = zeros(length(cmap_nums),3);
    for j=1:length(cmap_nums)
        time_cmap(ad_nums(j),:) = lineage_cmap(cmap_nums(j),:);
    end
    highlighted_time = create_highlighted_image(orig_i,ad_label_perim,'color_map',time_cmap);
    
    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
        edge_image_ad = create_highlighted_image(edge_image_ad,cell_edge,'color_map',edge_cmap(i_seen,:));
    end
    edge_image_ad = create_highlighted_image(edge_image_ad,bwperim(im2bw(ad_label,0)),'color_map',edge_cmap(i_seen,:));

    orig_i = orig_i(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3));
    highlighted_all = highlighted_all(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    highlighted_time = highlighted_time(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    edge_image_ad_bounded = edge_image_ad(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);

    spacer = 0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3);

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
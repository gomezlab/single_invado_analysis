function make_movie_frames(cfg_file)

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
i_size = size(imread(original_i_file,1));
i_count = size(imfinfo(original_i_file),2);
tracking_seqs = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect files to exclude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exclude_fullfile = fullfile(fileparts(original_i_file),excluded_frames_filename);
if (exist(exclude_fullfile,'file'))
    excluded_frames = load(exclude_fullfile);
else
    excluded_frames = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find extreme values in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extr_val_full_path = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_full_path,'file'))
    min_max = load(extr_val_full_path);
else
    min_max = find_extr_values(original_i_file);
    csvwrite(extr_val_full_path,min_max);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find extreme values in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bbox_full_path = fullfile(fileparts(original_i_file),bounding_box_filename);
if (exist(bbox_full_path,'file'))
    bounding_box = load(bbox_full_path);
else
    bounding_box = find_time_series_bbox(I_folder);
    csvwrite(bbox_full_path,bounding_box);
end

bounding_box(1:2) = bounding_box(1:2) - image_padding_min;
bounding_box(3:4) = bounding_box(3:4) + image_padding_min;
if (bounding_box(1) <= 0), bounding_box(1) = 1; end
if (bounding_box(2) <= 0), bounding_box(2) = 1; end
if (bounding_box(3) > i_size(2)), bounding_box(3) = i_size(2); end
if (bounding_box(4) > i_size(1)), bounding_box(4) = i_size(1); end


edge_c_map = jet(i_count);

pix_count = i_size(1)*i_size(2);
edge_image = ones(i_size(1),i_size(2),3);
edge_image_ad = ones(i_size(1),i_size(2),3);

max_live_adhesions = find_max_live_adhesions(tracking_seqs);

adhesion_tracking_map = jet(max_live_adhesions);

live_adhesion_to_color_map = zeros(size(tracking_seqs,1),1);

i_seen = 0;

for i = 1:i_count
    if (find(i==excluded_frames))
        continue;
    end
    if (i_seen + 1 > size(tracking_seqs,2))
        continue;
    end

    i_seen = i_seen + 1;
    
    orig_i = normalize_grayscale_image(imread(original_i_file,i),'min_max',min_max);

    padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
    padded_i_seen = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i_seen);

    adhesions = imread(fullfile(I_folder,padded_i_num,adhesions_filename));

    ad_label = bwlabel(adhesions,4);

    highlighted_all = cat(3,orig_i,orig_i,orig_i);
    search_ascending = 1;
    for j=1:size(tracking_seqs,1)
        if (tracking_seqs(j,i_seen) <= 0) 
            live_adhesion_to_color_map(j) = 0;
            continue;
        end
        
        if (live_adhesion_to_color_map(j))
            this_adhesion = zeros(i_size(1),i_size(2));
            this_adhesion(ad_label == tracking_seqs(j,i_seen)) = 1;
            this_color_map = adhesion_tracking_map(live_adhesion_to_color_map(j),:);
            highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_adhesion),'color',this_color_map);
        else
            short_list = live_adhesion_to_color_map(live_adhesion_to_color_map > 0);
            unique_used_colors = unique(short_list);
            
            if (search_ascending)
                search_seq = 1:max_live_adhesions;
                search_ascending = 0;
            else
                search_seq = max_live_adhesions:-1:1;
                search_ascending = 1;
            end
            
            poss_color = 0;
            for k = search_seq
                if (not(sum(unique_used_colors == k)))
                   poss_color = k;
                   break;
                end
            end
            
            if (not(poss_color))
                warning(['Could not find a unique color for adhesion lineage #: ', num2str(j)])
                live_adhesion_to_color_map(j) = 1;
            else        
                live_adhesion_to_color_map(j) = poss_color;
            end
            
            this_adhesion = zeros(i_size(1),i_size(2));
            this_adhesion(ad_label == tracking_seqs(j,i_seen)) = 1;
            this_color_map = adhesion_tracking_map(live_adhesion_to_color_map(j),:);
            highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_adhesion),'color',this_color_map);
        end
    end
    
    cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));

    for j = 1:3
        edge_image(find(cell_edge)+(j-1)*pix_count) = edge_c_map(i,j);
        edge_image_ad(find(cell_edge)+(j-1)*pix_count) = edge_c_map(i,j);
        edge_image_ad(find(bwperim(adhesions))+(j-1)*pix_count) = edge_c_map(i,j);
    end
    
    orig_i = orig_i(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3));
    highlighted_all = highlighted_all(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    edge_image_bounded = edge_image(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    edge_image_ad_bounded = edge_image_ad(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);

    frame = [cat(3,orig_i,orig_i,orig_i),0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3),highlighted_all];
    frame = {frame [edge_image_ad_bounded,0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3),highlighted_all]};
    
    if (exist('out_path','var'))
        for j = 1:length(out_prefix)
            output_filename = fullfile(out_path,[out_prefix{1,j},padded_i_seen,'.png']);
            fullpath = fileparts(output_filename);
            if (not(exist(fullpath,'dir')))
                mkdir(fullpath);
            end
            imwrite(frame{1,j},output_filename);
        end
    end
end

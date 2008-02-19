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
tracking_seqs = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect image numbers to exclude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist(excluded_frames_file,'file') && not(exist(excluded_frames_file,'dir')))
    excluded_frames = load(excluded_frames_file);
else
    excluded_frames = 0;
end

i_size = size(imread(fullfile(I_folder,sprintf(['%0',num2str(length(num2str(i_count))),'d'],1),focal_image)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist(bounding_box_file,'file'))
    bounding_box = load(bounding_box_file);
else
    bounding_box = find_time_series_bbox(I_folder);
    csvwrite(bbox_full_path,bounding_box_file);
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

    padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
    padded_i_seen = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i_seen);

    orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
    scale_factor = double(intmax(class(orig_i)));
    orig_i = double(orig_i)/scale_factor;
    
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
            
            if (i_seen > 1)
                highlighted_all = create_highlighted_image(highlighted_all,this_adhesion,'color',this_color_map);
            else
                highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_adhesion),'color',this_color_map);
            end
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
            output_filename = fullfile(out_path,out_prefix{1,j},[padded_i_seen,'.png']);
            fullpath = fileparts(output_filename);
            if (not(exist(fullpath,'dir')))
                mkdir(fullpath);
            end
            imwrite(frame{1,j},output_filename);
        end
    end
end

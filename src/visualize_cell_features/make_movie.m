function make_movie_frames(cfg_file)

[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

addpath(genpath(path_folders));

exclude_fullfile = fullfile(fileparts(original_i_file),excluded_frames_filename);
if (exist(exclude_fullfile,'file'))
    excluded_frames = load(exclude_fullfile);
else
    excluded_frames = 0;
end

if (exist('out_path','var'))
    if (not(exist(out_path,'dir')))
        mkdir(out_path);
    end
end

extr_val_fullfile = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_fullfile,'file'))
    min_max = load(extr_val_fullfile);
end

i_count = size(imfinfo(original_i_file),2);

edge_c_map = jet(i_count);

i_size = size(imread(original_i_file,1));
pix_count = i_size(1)*i_size(2);
edge_image = ones(i_size(1),i_size(2),3);
edge_image_ad = ones(i_size(1),i_size(2),3);

all_seqs = load(tracking_seq_file) + 1;

max_live_adhesions = find_max_live_adhesions(all_seqs);

adhesion_tracking_map = jet(max_live_adhesions);

live_adhesion_to_color_map = zeros(size(all_seqs,1));

i_seen = 0;

for i = 1:i_count
    if (find(i==excluded_frames))
        continue;
    end
    i_seen = i_seen + 1;
    
    orig_i = normalize_grayscale_image(imread(original_i_file,i),min_max(1),min_max(2));

    i_size = size(orig_i);

    padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
    padded_i_seen = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i_seen);

    I_1 = imread(fullfile(I_folder_1,padded_i_num,filename_1));
    
    I_2 = imread(fullfile(I_folder_2,padded_i_num,filename_2));
    markers = imread(fullfile(I_folder_1,padded_i_num,'focal_markers.png'));

    I_lab = bwlabel(I_2);
    I_test = zeros(i_size(1),i_size(2));
    I_test(find(I_lab == all_seqs(1,i_seen))) = 1;
    highlighted_test = create_highlighted_image(orig_i,bwperim(I_test));
    highlighted_test = create_highlighted_image(highlighted_test,bwperim(I_2 & not(I_test)),'color',1);    

    %highlighted_all = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_all = cat(3,orig_i,orig_i,orig_i);
    sort_ascending = 1;
    for j=1:size(all_seqs,1)
        if (not(all_seqs(j,i_seen))) 
            live_adhesion_to_color_map(j) = 0;
            continue;
        end
        
        if (live_adhesion_to_color_map(j))
            this_high = zeros(i_size(1),i_size(2));
            this_high(find(I_lab == all_seqs(j,i_seen))) = 1;
            this_color_map = adhesion_tracking_map(live_adhesion_to_color_map(j),:);
            highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_high),'color',this_color_map);
        else
            if (sort_ascending)
                sorted_used_colors = sort(unique(live_adhesion_to_color_map));
                sort_ascending = 0;
            else
                sorted_used_colors = sort(unique(live_adhesion_to_color_map),'descend');
                sort_ascending = 1;
            end
            
            poss_color = 0;
            for k = 1:max_live_adhesions
                find_l = length(find(sorted_used_colors == k));
                if (not(find_l))
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
            
            this_high = zeros(i_size(1),i_size(2));
            this_high(find(I_lab == all_seqs(j,i_seen))) = 1;
            this_color_map = adhesion_tracking_map(live_adhesion_to_color_map(j),:);
            highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_high),'color',this_color_map);
        end
        
    end
    
    highlighted_2 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(highlighted_1,bwperim(I_1 & not(I_2)),'color',1);
    %highlighted_1 = create_highlighted_image(highlighted_1,bwperim(markers),3);

    cell_edge_1 = bwperim(imread(fullfile(I_folder_1,padded_i_num,edge_filename)));

    for j = 1:3
        edge_image(find(cell_edge_1)+(j-1)*pix_count) = edge_c_map(i,j);
        edge_image_ad(find(cell_edge_1)+(j-1)*pix_count) = edge_c_map(i,j);
        edge_image_ad(find(bwperim(I_2))+(j-1)*pix_count) = edge_c_map(i,j);
    end

    
    %frame = [highlighted_1,0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_2];
    frame = [cat(3,orig_i,orig_i,orig_i),0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3),highlighted_all];
    frame = {frame [edge_image_ad,0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3),highlighted_all]};

    %frame = {[highlighted_1,0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3),highlighted_2]};
    
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

%system('/Users/mbergins/bin/ffmpeg -sameq -r 5 -b 8000 -y -i ../../results/norm/time_series_1/movie_files/time%03d.png -f huffyuv ../../results/norm/time_series_1/test.mp4');

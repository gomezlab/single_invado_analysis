function make_movie_frames(cfg_file,varargin)
profile on;
%MAKE_MOVIE_FRAMES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'MAKE_MOVIE_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});

debug = i_p.Results.debug;

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
tracking_seqs = load(tracking_seq_file) + 1;

i_size = size(imread(fullfile(I_folder,sprintf(['%0',num2str(length(num2str(i_count))),'d'],1),focal_image)));

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

edge_c_map = jet(i_count);

pix_count = i_size(1)*i_size(2);
edge_image_ad = ones(i_size(1),i_size(2),3);

max_live_adhesions = find_max_live_adhesions(tracking_seqs);

adhesion_tracking_cmap = jet(max_live_adhesions);
time_cmap = jet(size(tracking_seqs,2));

live_adhesion_to_color_map = zeros(size(tracking_seqs,1),1);
adhesion_time_to_cmap = zeros(size(tracking_seqs,1),1);

i_seen = 0;

for i = 1:i_count
    if (find(i==excluded_image_nums))
        continue;
    end
    if (i_seen + 1 > size(tracking_seqs,2))
        continue;
    end
    
    i_seen = i_seen + 1;
    
    if (debug && i_seen > 1); continue; end
    
    padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
    padded_i_seen = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i_seen);

    orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
    scale_factor = double(intmax(class(orig_i)));
    orig_i = double(orig_i)/scale_factor;
    
    adhesions = imread(fullfile(I_folder,padded_i_num,adhesions_filename));

    ad_label = bwlabel(adhesions,4);

    highlighted_all = cat(3,orig_i,orig_i,orig_i);
    
    %Uncomment the following line to highlight all the focal adhesions
    %with green before lineage colors are applied
    %highlighted_all = create_highlighted_image(highlighted_all,bwperim(adhesions));

    highlighted_time = highlighted_all;
    search_ascending = 1;
    
    any_highlights = 0;
    
    for j = 1:size(tracking_seqs,1)
        if (debug && j > 100); continue; end
        
        if (tracking_seqs(j,i_seen) <= 0)
            live_adhesion_to_color_map(j) = 0;
            continue;
        end
        
        any_highlights = 1;
        
        this_adhesion = zeros(i_size(1),i_size(2));
        this_adhesion(ad_label == tracking_seqs(j,i_seen)) = 1;
        
        %Time dependent colors
        if (adhesion_time_to_cmap(j))
            this_cmap = time_cmap(adhesion_time_to_cmap(j),:);
            highlighted_time = create_highlighted_image(highlighted_time,bwperim(this_adhesion),'color',this_cmap);
        else
            adhesion_time_to_cmap(j) = i_seen;
            this_cmap = time_cmap(adhesion_time_to_cmap(j),:);
            highlighted_time = create_highlighted_image(highlighted_time,bwperim(this_adhesion),'color',this_cmap);
        end

        %Unique adhesion colors
        if (live_adhesion_to_color_map(j))
            this_cmap = adhesion_tracking_cmap(live_adhesion_to_color_map(j),:);
            highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_adhesion),'color',this_cmap);
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
                warning('Matlab:adColor',['Could not find a unique color for adhesion lineage #: ', num2str(j)])
                live_adhesion_to_color_map(j) = 1;
            else        
                live_adhesion_to_color_map(j) = poss_color;
            end
            
            this_cmap = adhesion_tracking_cmap(live_adhesion_to_color_map(j),:);
            
            if (i_seen > 1)
                highlighted_all = create_highlighted_image(highlighted_all,this_adhesion,'color',this_cmap);
            else
                highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_adhesion),'color',this_cmap);
            end
        end
    end

    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
        for j = 1:3
            edge_image_ad(find(cell_edge)+(j-1)*pix_count) = edge_c_map(i,j);
        end
    else
        for j = 1:3
            edge_image_ad(find(bwperim(adhesions))+(j-1)*pix_count) = edge_c_map(i,j);
        end
    end
    
    orig_i = orig_i(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3));
    highlighted_all = highlighted_all(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    highlighted_time = highlighted_time(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    edge_image_ad_bounded = edge_image_ad(bounding_box(2):bounding_box(4), bounding_box(1):bounding_box(3),1:3);
    
    spacer = 0.5*ones(size(orig_i,1),round(0.02*size(orig_i,2)),3);
    
    if (not(any_highlights)); continue; end
    
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
end

profile off

if (debug)
    profile viewer
end
function make_movie_frames(cfg_file,varargin)
%MAKE_MOVIE_FRAMES    Builds movie frames with the adhesions highlighted in
%                     using various conventions
%
%   make_movie_frames(cfg_file,options) builds individual movie frames from
%   raw experimental data, where files are placed and the movie config
%   options are coded in cfg_file
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'MAKE_MOVIE_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('no_scale_bar',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(cfg_file);
while 1
    line = fgetl(fid);
    if ~ischar(line), break; end
    eval(line);
end

addpath(genpath(path_folders));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect General Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenames = add_filenames_to_struct(struct());

max_image_num = find_max_image_num(I_folder);
folder_char_length = length(num2str(max_image_num));
i_size = size(imread(fullfile(I_folder,num2str(max_image_num),filenames.puncta_filename)));

tracking_seq = load(tracking_seq_file) + 1;
%occastionally, there will be fields were no cells were detected, in that
%case, the tracking matrix is an empty file, but we still need a matrix to
%index into, so make an empty tracking matrix with a line of zeros
if (isempty(tracking_seq))
    tracking_seq = zeros(1,max_image_num);
end

max_live_adhesions = find_max_live_adhesions(tracking_seq);
lineage_cmap = jet(max_live_adhesions);

lineage_to_cmap_row = assign_unique_colors_to_lineages(tracking_seq);

for i = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather and scale the input puncta image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_image = double(imread(fullfile(I_folder,padded_i_num,filenames.puncta_filename)));
    puncta_limits = csvread(fullfile(I_folder,padded_i_num,filenames.puncta_range_file));
    
    puncta_image = puncta_image - puncta_limits(1);
    puncta_image = puncta_image .* (1/puncta_limits(2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather and scale the gel image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_image = double(imread(fullfile(I_folder,padded_i_num,filenames.gel_filename)));
    gel_limits = csvread(fullfile(I_folder,padded_i_num,filenames.gel_range_file));
    
    gel_image = gel_image - gel_limits(1);
    gel_image = gel_image .* (1/gel_limits(2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_mask_label = imread(fullfile(I_folder,padded_i_num,filenames.labeled_cell_mask_filename));
    cell_mask_label_perim = zeros(size(cell_mask_label));
    for j=1:max(cell_mask_label(:))
        this_cell = zeros(size(cell_mask_label));
        this_cell(cell_mask_label == j) = 1;
        this_cell = bwperim(this_cell);
        cell_mask_label_perim(this_cell) = j;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    live_rows = tracking_seq(:,i) > 0;
    ad_nums = tracking_seq(live_rows,i);
    
    %Build the unique lineage highlighted image
    this_cmap(ad_nums,:) = lineage_cmap(lineage_to_cmap_row(live_rows),:); %#ok<AGROW>
    highlighted_puncta = create_highlighted_image(puncta_image,cell_mask_label_perim,'color_map',this_cmap);
    
%     highlighted_gel = create_highlighted_image(gel_image,cell_mask_label_perim,'color_map',this_cmap);
    
    if (not(exist(out_path,'dir')))
        mkdir(out_path);
    end
    
%     imwrite(highlighted_puncta,fullfile(out_path,[padded_i_num,'.png']));
    imwrite(highlighted_puncta,fullfile(I_folder,padded_i_num,'tracking.png'));
%
%     %Build the birth time highlighted image
%     cmap_nums = birth_time_to_cmap(tracking_seq(:,i_seen) > 0);
%     assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in birth time image %d',padded_i_num);
%     this_cmap = zeros(max(puncta_label_perim(:)),3);
%     this_cmap(ad_nums,:) = time_cmap(cmap_nums,:);
%     highlighted_time = create_highlighted_image(puncta_image,puncta_label_perim,'color_map',this_cmap);
%
%     edge_image_ad = create_highlighted_image(edge_image_ad,im2bw(puncta_label_perim,0),'color_map',edge_cmap(i_seen,:));
%
%     spacer = 0.5*ones(size(edge_image_ad_bounded,1),1,3);
%
%     frame = cell(1,3);
%     frame{1} = [edge_image_ad_bounded,spacer,highlighted_all];
%     frame{2} = [cat(3,puncta_image,puncta_image,puncta_image),spacer,highlighted_all];
%     frame{3} = [edge_image_ad_bounded,spacer,highlighted_time];
%
%     %Add scale bars if the pixel size is available
%     if (exist('pixel_size','var') && not(i_p.Results.no_scale_bar))
%         for j = 1:size(frame,2)
%             frame{j} = draw_scale_bar(frame{j},pixel_size);
%         end
%     end
end

toc;

function lineage_cmap = assign_unique_colors_to_lineages(tracking_seq)

max_live_adhesions = find_max_live_adhesions(tracking_seq);
lineage_cmap = zeros(size(tracking_seq,1),1);

use_beginning = 0;

for time = 1:size(tracking_seq,2)
    living_lineage_rows = find(tracking_seq(:,time) > 0);
    
    for j = 1:length(living_lineage_rows)        
        this_row = living_lineage_rows(j);
        
        %if a color number is already assigned, skip to the next row with a
        %living cell
        if (lineage_cmap(this_row) > 0)
            continue;
        end
        
        %find all the currently in use color numbers
        used_c_nums = sort(lineage_cmap(living_lineage_rows));
        
        %throw out all the zeros
        used_c_nums = used_c_nums(used_c_nums ~= 0);
        
        unused_c_nums = setdiff(1:max_live_adhesions,used_c_nums);
        
        assert(not(isempty(unused_c_nums)), 'Couldn''t find an unused color number')
        
        if (use_beginning == 0)
            lineage_cmap(this_row) = unused_c_nums(1);
            use_beginning = 1;
        else
            lineage_cmap(this_row) = unused_c_nums(end);
            use_beginning = 0;
        end
        1;
    end
end

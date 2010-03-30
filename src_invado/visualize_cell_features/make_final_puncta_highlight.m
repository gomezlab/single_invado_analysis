function make_final_puncta_highlight(cfg_file,invado_file,areas_lineage_ts_file,varargin)
%MAKE_FINAL_PUNCTA_HIGHLIGHT    Builds several still image visualizations
%                               of the invadopodia in a given experiment,
%                               highlighting them over the final image of
%                               the visualization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'MAKE_FINAL_PUNCTA_HIGHLIGHT';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addRequired('invado_file',@(x)exist(x,'file') == 2);
i_p.addRequired('areas_lineage_ts_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('no_scale_bar',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,invado_file,areas_lineage_ts_file,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

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
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invado_data = csvread(invado_file);

%after loading the tracking sequence filter to only include those puncta
%included in the invadopodia list, remember the list formated so the first
%column contains the lineage number, with the first lineage as 1, so no
%need to translate
tracking_seq = load(tracking_seq_file) + 1;
tracking_seq = tracking_seq(invado_data(:,1),:);

%read in the area lineage file and filter down to only the included
%lineages
areas_lineage_ts = csvread(areas_lineage_ts_file);
areas_lineage_ts = areas_lineage_ts(invado_data(:,1),:);

for i = 1:size(areas_lineage_ts,1)
    largest_area_index = find(max(areas_lineage_ts(i,:)) == areas_lineage_ts(i,:),1);
    temp_row = zeros(1,size(tracking_seq,2));
    temp_row(largest_area_index) = tracking_seq(i,largest_area_index);
    tracking_seq(i,:) = temp_row;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge_cmap = jet(size(tracking_seq,2));
% %define the edge image here because the edge image will be added to each
% %image loop, so the image should be global
% edge_image_ad = ones(i_size(1),i_size(2),3);

lineage_cmap = jet(size(tracking_seq,1));
lineage_to_cmap = 1:size(tracking_seq,1);

% time_cmap = hot(size(tracking_seq,2));
% birth_time_to_cmap = zeros(size(tracking_seq,1),1);

%fill out the p_value based color map
p_value_ranges = linspace(0.005,0.05,15);
p_value_cmap = hot(length(p_value_ranges));
%flip the rows of the color map
p_value_cmap = p_value_cmap(size(p_value_cmap,1):-1:1,:);
lineage_to_p_value_cmap = zeros(size(tracking_seq,1),1);
for i = 1:length(invado_data(:,2))
    try
        lineage_to_p_value_cmap(i) = find(p_value_ranges >= invado_data(i,2), 1);
    catch
        lineage_to_p_value_cmap(i) = length(p_value_ranges);
    end
end

%fill out the local_diff based color map
local_diff_ranges = linspace(-0.01,0,15);
local_diff_cmap = cool(length(local_diff_ranges));
%flip the rows of the color map
% local_diff_cmap = local_diff_cmap(size(local_diff_cmap,1):-1:1,:);
lineage_to_local_diff_cmap = zeros(size(tracking_seq,1),1);
for i = 1:length(invado_data(:,3))
    try
        lineage_to_local_diff_cmap(i) = find(local_diff_ranges >= invado_data(i,3), 1);
    catch
        lineage_to_local_diff_cmap(i) = length(local_diff_ranges);
    end
end

%final all the image directories
image_dirs = dir(I_folder);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%read in the last gel image and normalize to 0-1
last_gel_image = double(imread(fullfile(I_folder,image_dirs(end).name,gel_image)));
last_binary_shift = logical(imread(fullfile(I_folder,image_dirs(end).name,'binary_shift.png')));
gel_limits = csvread(fullfile(I_folder,image_dirs(end).name,'gel_image_range.csv'));
last_gel_image = last_gel_image - gel_limits(1);
last_gel_image = last_gel_image .* (1/gel_limits(2));
last_gel_image(not(last_binary_shift)) = 0;
last_gel_image = cat(3,last_gel_image,last_gel_image,last_gel_image);

last_gel_image_high = last_gel_image;
last_gel_image_p_val_high = last_gel_image;
last_gel_image_local_diff_high = last_gel_image;

for i = 1:length(image_dirs)
    %check for the presence of adhesions to map onto the last gel image, if
    %not present, skip to the next image folder
    if (not(any(tracking_seq(:,i) > 0)))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label_perim = imread(fullfile(I_folder,image_dirs(i).name,adhesions_perim_filename));
    
    ad_nums = tracking_seq(:,i);
    ad_nums = ad_nums(ad_nums > 0);
    
    temp_ad_label_perim = zeros(size(ad_label_perim));
    temp_ad_label_perim(ismember(ad_label_perim,ad_nums)) = ad_label_perim(ismember(ad_label_perim,ad_nums));
    
    ad_label_perim = temp_ad_label_perim;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',i);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = lineage_cmap(cmap_nums,:);
    last_gel_image_high = create_highlighted_image(last_gel_image_high,ad_label_perim,'color_map',this_cmap, 'mix_percent',0.5);
    
    %Build the p-value colored highlight image
    cmap_nums = lineage_to_p_value_cmap(tracking_seq(:,i) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',i);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = p_value_cmap(cmap_nums,:);
    last_gel_image_p_val_high = create_highlighted_image(last_gel_image_p_val_high,ad_label_perim,'color_map',this_cmap, 'mix_percent',0.5);

    %Build the local difference colored highlight image
    cmap_nums = lineage_to_local_diff_cmap(tracking_seq(:,i) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',i);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = local_diff_cmap(cmap_nums,:);
    last_gel_image_local_diff_high = create_highlighted_image(last_gel_image_local_diff_high,ad_label_perim,'color_map',this_cmap, 'mix_percent',0.5);    
    
    if(i_p.Results.debug), disp(i); end
end

%build and output the images
spacer = ones(size(last_gel_image,1),1,3);

if (exist('pixel_size','var') && not(i_p.Results.no_scale_bar))
    last_gel_image_high = draw_scale_bar(last_gel_image_high,pixel_size);
    last_gel_image_p_val_high = draw_scale_bar(last_gel_image_p_val_high,pixel_size);
end

imwrite([last_gel_image, spacer, last_gel_image_high], fullfile(vis_folder,'invado_highlight_last.png'))
imwrite([last_gel_image, spacer, last_gel_image_p_val_high], fullfile(vis_folder,'invado_highlight_p_val.png'))
imwrite([last_gel_image, spacer, last_gel_image_local_diff_high], fullfile(vis_folder,'invado_highlight_local_diff.png'))

spacer_long = ones(1,size(last_gel_image,2)*2+1,3);

imwrite([last_gel_image, spacer, last_gel_image_high;spacer_long;last_gel_image_p_val_high,spacer,last_gel_image_local_diff_high], fullfile(vis_folder,'invado_highlight_all.png'))

profile off;
if (i_p.Results.debug), profile viewer; end

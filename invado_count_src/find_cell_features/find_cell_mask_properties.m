function find_cell_mask_properties(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

no_cell_diffs = [];
overlap_diffs = [];
data_points = 0;
overlap_sizes = [];
cell_diffs = [];
p_vals = [];

for i_num = 1:size(image_dirs,1)
    if (i_num > 1)
        prior_dir = fullfile(base_dir,image_dirs(i_num-1).name);
        if (exist('current_data','var'))
            prior_data = current_data;
        else
            prior_data = read_in_file_set(prior_dir,filenames);
        end
    else
        prior_dir = fullfile(base_dir,image_dirs(i_num).name);
        prior_data = read_in_file_set(prior_dir,filenames);        
    end
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    cell_props = collect_cell_properties(current_data,prior_data,'debug',i_p.Results.debug);
    
    %convert cell props data into movieInfo data suitable for u-track
%     [xCoord, yCoord, amp] = convert_props_to_movieInfo(cell_props);
%     movieInfo(i_num).xCoord = xCoord;
%     movieInfo(i_num).yCoord = yCoord;
%     movieInfo(i_num).amp = amp;

    data_points = data_points + length(cell_props);
    tracking_props = collect_tracking_properties(current_data,prior_data,'debug',i_p.Results.debug);
    if (i_num ~= 1)
        overlap_sizes = [overlap_sizes, [cell_props.Overlap_region_size]];
        cell_diffs = [cell_diffs, [cell_props.Cell_gel_diff]];
        no_cell_diffs = [no_cell_diffs, [cell_props(1).no_cells_diff]'];
        p_vals = [p_vals, [cell_props.Cell_gel_diff_p_val]];
        
        cell_props = rmfield(cell_props,'no_cells_diff');
    end

    write_adhesion_data(cell_props,'out_dir',fullfile(current_dir,'raw_data'));
    write_adhesion_data(tracking_props,'out_dir',fullfile(prior_dir,'raw_data'));
    
    if (i_num == size(image_dirs,1))
        tracking_props = collect_tracking_properties(current_data,current_data,'debug',i_p.Results.debug);
        write_adhesion_data(tracking_props,'out_dir',fullfile(current_dir,'raw_data'));
    end
    
    %make a diagnostic figure showing where the cell masks overlap with the
    %next image
    temp = double(current_data.cell_mask);
    temp(prior_data.cell_mask) = 2;
    temp(prior_data.cell_mask & current_data.cell_mask) = 3;
    imwrite(label2rgb(temp),fullfile(current_dir,'cell_overlaps.png'));

    %make an image showing where the cells are located overlayed on the gel
    %image
    gel_range_norm = create_highlighted_image(current_data.gel_image_norm,bwperim(current_data.cell_mask));
    imwrite(gel_range_norm,fullfile(current_dir,'gel_highlights.png'));

    if (mod(i_num,10)==0)
        disp(['Done processing image number: ',num2str(i_num)])
    end    
end

% sample_sizes = 1000:1000:20000;
% min_maxes = NaN(2,length(sample_sizes));
% for i=1:length(sample_sizes)
%     min_maxes(:,i) = bootstrap_mean_small(no_cell_diffs,'samp_size',sample_sizes(i))';
%     disp(sample_sizes(i));
%     plot(sample_sizes, min_maxes);
% end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xCoord,yCoord,amp] = convert_props_to_movieInfo(cell_props)

xCoord = zeros(length(cell_props),2);
yCoord = zeros(length(cell_props),2);
amp = zeros(length(cell_props),2);

for i=1:length(cell_props)
    xCoord(i,1) = cell_props(i).Centroid(1);
    yCoord(i,1) = cell_props(i).Centroid(2);
    
    amp(i,1) = cell_props(i).MeanIntensity;
    amp(i,2) = cell_props(i).StdIntensity;
end

function cell_props = collect_cell_properties(current_data,prior_data,varargin)
% COLLECT_cell_PROPERTIES    using the identified cells, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('current_data',@isstruct);
i_p.addRequired('prior_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prior_data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_props = regionprops(current_data.labeled_cells,current_data.gel_image, ...
    'Area','Centroid','MeanIntensity','PixelValues');

for i=1:max(current_data.labeled_cells(:))
        cell_props(i).StdIntensity = std([cell_props(i).PixelValues]);
end
cell_props = rmfield(cell_props,'PixelValues');

[cell_props.Overlap_region_size] = deal(NaN);
[cell_props.Cell_gel_diff_p_val] = deal(NaN);
[cell_props.Cell_gel_diff] = deal(NaN);
[cell_props.Cell_gel_diff_median] = deal(NaN);
[cell_props.Cell_gel_diff_total] = deal(NaN);

%when the first image is both the prior and current data, we only want the
%area and centroid variables to contain any data, the other properties are
%not relavent for the first image
if (all(current_data.gel_image(:) == prior_data.gel_image(:)) == 1)
    return;
end

cell_props(1).no_cells_diff = current_data.gel_image(current_data.no_cells)*current_data.intensity_correction - ...
    prior_data.gel_image(prior_data.no_cells)*prior_data.intensity_correction;
% cell_props(1).no_cells_diff = current_data.gel_image(current_data.no_cells) - ...
%     prior_data.gel_image(prior_data.no_cells);
% [h,p,ci] = ttest(cell_props(1).no_cells_diff);
% ci

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(current_data.labeled_cells(:))
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_cell = current_data.labeled_cells;
    this_cell(current_data.labeled_cells ~= i) = 0;
    this_cell = logical(this_cell);
    
    prev_cells = prior_data.cell_mask;
    
    overlap_region = this_cell & prev_cells;
    
    differences = current_data.gel_image(overlap_region)*current_data.intensity_correction - ... 
        prior_data.gel_image(overlap_region)*prior_data.intensity_correction;    
    
    [h,p] = ttest(differences);
    cell_props(i).Cell_gel_diff = mean(differences);
    cell_props(i).Cell_gel_diff_p_val = p;
    cell_props(i).Cell_gel_diff_median = median(differences);
    cell_props(i).Cell_gel_diff_total = sum(differences);
    
    %single cell diagnostics
    if (i_p.Results.debug)
        temp = double(this_cell);
        temp(prev_cells) = 2;
        temp(overlap_region) = 3;
        subplot(1,2,1); imshow(label2rgb(temp));
        subplot(1,2,2); hist(differences);
        1;
%         imshow(current_data.gel_image.*overlap_region - prior_data.gel_image.*overlap_region,[]);
    end
end

function prior_props = collect_tracking_properties(current_data,prior_data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('current_data',@isstruct);
i_p.addRequired('prior_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prior_data,varargin{:});

current_props = regionprops(current_data.labeled_cells,'Centroid');
prior_props = regionprops(prior_data.labeled_cells,'Centroid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(prior_data.labeled_cells(:))
    for j=1:max(current_data.labeled_cells(:))
        prior_props(i).Cent_dist(j) = sqrt((current_props(j).Centroid(1) - prior_props(i).Centroid(1))^2 + ...
            (current_props(j).Centroid(2) - prior_props(i).Centroid(2))^2);
    end
end

for i=1:max(prior_data.labeled_cells(:))
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_cell = prior_data.labeled_cells;
    this_cell(prior_data.labeled_cells ~= i) = 0;
    this_cell = logical(this_cell);
    
    prior_props(i).Pix_sim = zeros(max(current_data.labeled_cells(:)),1);
    this_cell_area = sum(sum(this_cell));
    
    next_overlaping_cells = current_data.labeled_cells(this_cell);
    
    overlap_nums = unique(next_overlaping_cells);
    if (overlap_nums(1) == 0)
        overlap_nums = overlap_nums(2:end);
    end
    for j=1:length(overlap_nums)
       this_overlap_num = overlap_nums(j);
       
       overlap_area = sum(next_overlaping_cells == this_overlap_num);
       prior_props(i).Pix_sim(this_overlap_num) = overlap_area/this_cell_area;
    end
end

%we don't want to write the centroid field back out, so remove it from the
%structure
prior_props = rmfield(prior_props,'Centroid');


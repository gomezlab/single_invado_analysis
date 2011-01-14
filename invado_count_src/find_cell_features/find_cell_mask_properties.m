function find_cell_mask_properties(exp_dir,varargin)

tic; profile off; profile on;
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
    tracking_props = collect_tracking_properties(current_data,prior_data,'debug',i_p.Results.debug);

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
    gel_range_norm = double(current_data.gel_no_norm) - current_data.gel_range(1);
    gel_range_norm = gel_range_norm / (current_data.gel_range(2) - current_data.gel_range(1));
    gel_range_norm = create_highlighted_image(gel_range_norm,bwperim(current_data.cell_mask));
    imwrite(gel_range_norm,fullfile(current_dir,'gel_highlights.png'));

    if (mod(i_num,10)==0)
        disp(['Done processing image number: ',num2str(i_num)])
    end    
end

toc;
profile off;
if (i_p.Results.debug), profile viewer; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cell_props = collect_cell_properties(current_data,prior_data,varargin)
% COLLECT_cell_PROPERTIES    using the identified cells, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_CELL_PROPERTIES';

i_p.addRequired('current_data',@isstruct);
i_p.addRequired('prior_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prior_data,varargin{:});

cell_props = regionprops(current_data.labeled_cells,'Area','Centroid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    cell_props(i).Overlap_region_size = sum(sum(overlap_region));
    
    differences = current_data.gel_no_norm(overlap_region)*current_data.intensity_correction - ... 
        prior_data.gel_no_norm(overlap_region)*prior_data.intensity_correction;
    
    [h,p] = ttest(differences);
    cell_props(i).Cell_gel_diff_p_val = p;
    cell_props(i).Cell_gel_diff = mean(differences);
    cell_props(i).Cell_gel_diff_percent = mean(differences)/max(current_data.gel_range);
    
    %single cell diagnostics
    if (i_p.Results.debug)
        temp = double(this_cell);
        temp(prev_cells) = 2;
        temp(overlap_region) = 3;
        subplot(1,2,1); imshow(label2rgb(temp));
        subplot(1,2,2); hist(differences);
        1;
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


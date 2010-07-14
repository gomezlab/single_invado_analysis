function find_cell_mask_properties(current_dir,prev_dir,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;

i_p = inputParser;

i_p.addRequired('current_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('prev_dir',@(x)exist(x,'dir') == 7);

i_p.parse(current_dir,prev_dir);

i_p.addParamValue('output_dir',current_dir,@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(current_dir,prev_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_data = struct;

%read in and normalize the input gel image
current_data.gel_image  = imread(fullfile(current_dir, filenames.gel_filename));
current_data.gel_no_norm = current_data.gel_image;
scale_factor = double(intmax(class(current_data.gel_image)));
current_data.gel_image  = double(current_data.gel_image)/scale_factor;

%read in gel min/max
current_data.gel_range = csvread(fullfile(current_dir, filenames.gel_range_file));

%read in the intensity correction coefficient
current_data.intensity_correction = csvread(fullfile(current_dir, filenames.intensity_correction_filename));

%read in the cell mask file
current_data.cell_mask = logical(imread(fullfile(current_dir, filenames.cell_mask_filename)));

%read in the labeled cells
current_data.labeled_cells = imread(fullfile(current_dir, filenames.labeled_cell_mask_filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the previous directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prev_data = struct;

prev_data.gel_image  = imread(fullfile(prev_dir, filenames.gel_filename));
scale_factor = double(intmax(class(prev_data.gel_image)));
prev_data.gel_image  = double(prev_data.gel_image)/scale_factor;

prev_data.binary_shift = logical(imread(fullfile(prev_dir, filenames.binary_shift_filename)));
prev_data.intensity_correction = csvread(fullfile(prev_dir, filenames.intensity_correction_filename));
prev_data.cell_mask = logical(imread(fullfile(prev_dir, filenames.cell_mask_filename)));
prev_data.labeled_cells = imread(fullfile(prev_dir, filenames.labeled_cell_mask_filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = double(current_data.cell_mask);
temp(prev_data.cell_mask) = 2;
temp(prev_data.cell_mask & current_data.cell_mask) = 3;
imwrite(label2rgb(temp),fullfile(i_p.Results.output_dir,'cell_overlaps.png'));



cell_properties = collect_cell_properties(current_data,prev_data,'debug',i_p.Results.debug);

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(cell_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));

profile off;
if (i_p.Results.debug), profile viewer; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cell_props = collect_cell_properties(current_data,prev_data,varargin)
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
i_p.addRequired('prev_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prev_data,varargin{:});

cell_props = regionprops(current_data.labeled_cells,'all');

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
    
    prev_cells = prev_data.cell_mask;
    
    overlap_region = this_cell & prev_cells;
        
    differences = current_data.gel_image(overlap_region)*current_data.intensity_correction - prev_data.gel_image(overlap_region)*prev_data.intensity_correction;
%     differences = current_data.gel_image(this_cell) - prev_data.gel_image(this_cell);
    [h,p] = ttest(differences);
    
    cell_props(i).Cell_gel_diff_p_val = p;
    cell_props(i).Cell_gel_diff = mean(differences);
end
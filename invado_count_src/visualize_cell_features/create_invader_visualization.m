function create_invader_visualization(current_dir,first_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;

i_p = inputParser;

i_p.addRequired('first_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('current_dir',@(x)exist(x,'dir') == 7);

i_p.parse(current_dir, first_dir);

i_p.addParamValue('output_dir',current_dir,@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(current_dir, first_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('../find_cell_features/matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_data = struct;

current_data.binary_shift  = logical(imread(fullfile(current_dir, filenames.binary_shift_filename)));

%read in and normalize the input focal adhesion image
current_data.puncta_image  = imread(fullfile(current_dir, filenames.puncta_filename));
current_data.puncta_range  = csvread(fullfile(current_dir, filenames.puncta_range_file));
current_data.puncta_image = double(current_data.puncta_image) - current_data.puncta_range(1);
current_data.puncta_image = current_data.puncta_image ./ (current_data.puncta_range(2) - current_data.puncta_range(1));
current_data.puncta_image = current_data.puncta_image .* current_data.binary_shift;

%read in and normalize the input focal adhesion image
current_data.gel_image  = imread(fullfile(current_dir, filenames.gel_filename));
current_data.gel_range  = csvread(fullfile(current_dir, filenames.gel_range_file));
current_data.gel_image = double(current_data.gel_image) - current_data.gel_range(1);
current_data.gel_image = current_data.gel_image ./ (current_data.gel_range(2) - current_data.gel_range(1));
current_data.gel_image = current_data.gel_image .* current_data.binary_shift;

%read in the intensity correction coefficient
current_data.intensity_correction = csvread(fullfile(current_dir, filenames.intensity_correction_filename));

%read in the cell mask file
current_data.cell_mask = logical(imread(fullfile(current_dir, filenames.cell_mask_filename)));

%read in the labeled cells
current_data.labeled_cell_mask = imread(fullfile(current_dir, filenames.labeled_cell_mask_filename));

current_data.Cell_diff = csvread(fullfile(current_dir, 'raw_data','Cell_gel_diff.csv'));
current_data.Cell_diff_p_val = csvread(fullfile(current_dir, 'raw_data','Cell_gel_diff_p_val.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the first directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_data = struct;

%read in and normalize the input focal adhesion image
first_data.gel_image  = imread(fullfile(first_dir, filenames.gel_filename));
first_data.gel_image = double(first_data.gel_image) - current_data.gel_range(1);
first_data.gel_image = first_data.gel_image ./ (current_data.gel_range(2) - current_data.gel_range(1));

invader_thresh = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine the number of comparisons made in this experimental set and
%correct the invader threshold by the Bonferroni correction
comparison_number = determine_comparison_number(fullfile(current_dir,'..','..','..'));
invader_thresh = invader_thresh/comparison_number;

non_invaders = ismember(current_data.labeled_cell_mask, find(not(isnan(current_data.Cell_diff_p_val)) & (current_data.Cell_diff_p_val > invader_thresh | current_data.Cell_diff > 0)));
invaders = ismember(current_data.labeled_cell_mask, find(current_data.Cell_diff < 0 & current_data.Cell_diff_p_val < invader_thresh));
no_class = ismember(current_data.labeled_cell_mask, find(isnan(current_data.Cell_diff_p_val)));

assert(sum(sum(non_invaders & invaders)) == 0)
assert(sum(sum(non_invaders & no_class)) == 0)
assert(sum(sum(invaders & no_class)) == 0)

assert(sum(sum(non_invaders | no_class | invaders)) == sum(sum(current_data.cell_mask)))

%we only need to clear out the area outside the binary shift of the current
%data image because the first directory binary shift should also include
%the entire image
diff_image = (current_data.gel_image - first_data.gel_image).*current_data.binary_shift;
diff_image = diff_image - min(diff_image(:));
diff_image = diff_image ./ (max(diff_image(:)) - min(diff_image(:)));

diff_image = create_highlighted_image(diff_image,bwperim(no_class),'color_map',[0,0,1]);
diff_image = create_highlighted_image(diff_image,bwperim(non_invaders),'color_map',[1,0,0]);
diff_image = create_highlighted_image(diff_image,bwperim(invaders),'color_map',[0,1,0]);

imwrite(diff_image, fullfile(current_dir, 'diff_invader.png'))

gel_highlight = create_highlighted_image(current_data.gel_image,bwperim(no_class),'color_map',[0,0,1]);
gel_highlight = create_highlighted_image(gel_highlight,bwperim(non_invaders),'color_map',[1,0,0]);
gel_highlight = create_highlighted_image(gel_highlight,bwperim(invaders),'color_map',[0,1,0]);

puncta_highlight = create_highlighted_image(current_data.puncta_image,bwperim(no_class),'color_map',[0,0,1]);
puncta_highlight = create_highlighted_image(puncta_highlight,bwperim(non_invaders),'color_map',[1,0,0]);
puncta_highlight = create_highlighted_image(puncta_highlight,bwperim(invaders),'color_map',[0,1,0]);

spacer = ones(size(puncta_highlight,1),1,3);

% imwrite([gel_highlight,spacer,puncta_highlight, spacer, diff_image], fullfile(current_dir, 'invader_and_not.png'))
imwrite([gel_highlight,spacer,puncta_highlight], fullfile(current_dir, 'invader_and_not.png'))

function comp_number = determine_comparison_number(experiment_dir)

single_exp_dirs = dir(experiment_dir);

assert(single_exp_dirs(1).name == '.')
assert(all(single_exp_dirs(2).name=='..'))
single_exp_dirs = single_exp_dirs(3:end);

p_vals = [];

for i=1:length(single_exp_dirs)
    data_file = fullfile(experiment_dir,single_exp_dirs(i).name, ...
        'adhesion_props','lin_time_series','Cell_gel_diff_p_val.csv');
    if (not(exist(data_file,'file')))
        warn(['Can''t find: ',data_file])
    end
    
    load(data_file);
    p_vals = [p_vals; Cell_gel_diff_p_val];
end

comp_number = sum(sum(not(isnan(p_vals))));
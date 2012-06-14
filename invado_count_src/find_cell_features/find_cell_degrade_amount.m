function find_cell_degrade_amount(exp_dir,varargin)

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

tracking_file = fullfile(base_dir, image_dirs(1).name,filenames.tracking);
%check for the existance of a tracking file, if absent, there weren't any
%cells in this field, return from the function
if (not(exist(tracking_file,'file')))
    disp('No tracking matrix found, assuming no cells in field');
    return;
else
    tracking_mat = csvread(tracking_file)+1;
end

final_data = read_in_file_set(fullfile(base_dir,image_dirs(size(image_dirs,1)).name),filenames);
first_data = read_in_file_set(fullfile(base_dir,image_dirs(1).name),filenames);

cell_hit_counts =  cell(size(tracking_mat,1),1);

for i_num = 1:size(tracking_mat,2)
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    for cell_num = 1:size(tracking_mat,1)
        if (tracking_mat(cell_num,i_num) == 0), continue; end
        
        if (isempty(cell_hit_counts{cell_num})) 
            cell_hit_counts{cell_num} = zeros(size(current_data.labeled_cells,1),size(current_data.labeled_cells,2));
        end
        
        cell_hit_counts{cell_num} = cell_hit_counts{cell_num} + double(current_data.labeled_cells == tracking_mat(cell_num,i_num));
    end
end

diffs = zeros(length(cell_hit_counts),1);

for cell_num = 1:length(cell_hit_counts)
    cell_extent = cell_hit_counts{cell_num} > 10;
    diffs(cell_num) = mean(final_data.gel_image(cell_extent) - first_data.gel_image(cell_extent))/386;
end

csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.final_gel_diffs),diffs);

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
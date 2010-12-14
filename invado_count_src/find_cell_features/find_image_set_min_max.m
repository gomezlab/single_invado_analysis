function find_image_set_min_max(base_dir, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_IMAGE_SET_MIN_MAX';

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

addpath('matlab_scripts');

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

gel_image_range = [Inf -Inf];
puncta_image_range = [Inf -Inf];

for i = 1:size(image_dirs)
    gel_image = imread(fullfile(base_dir,image_dirs(i).name,filenames.gel_filename));
    this_gel_min_max = [min(gel_image(:)), max(gel_image(:))];
    if (this_gel_min_max(1) < gel_image_range(1))
        gel_image_range(1) = this_gel_min_max(1);
    end
    if (this_gel_min_max(2) > gel_image_range(2))
        gel_image_range(2) = this_gel_min_max(2);
    end
        
    puncta_image = imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta_filename));
    this_puncta_min_max = [min(puncta_image(:)), max(puncta_image(:))];
    if (this_puncta_min_max(1) < puncta_image_range(1))
        puncta_image_range(1) = this_puncta_min_max(1);
    end
    if (this_puncta_min_max(2) > puncta_image_range(2))
        puncta_image_range(2) = this_puncta_min_max(2);
    end
end

for i = 1:size(image_dirs)
    csvwrite(fullfile(base_dir,image_dirs(i).name,filenames.gel_range_file),gel_image_range)
    csvwrite(fullfile(base_dir,image_dirs(i).name,filenames.puncta_range_file),puncta_image_range)
end

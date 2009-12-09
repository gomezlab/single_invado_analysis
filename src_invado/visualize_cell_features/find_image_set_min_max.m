function find_image_set_min_max(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_IMAGE_SET_MIN_MAX';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_folder,varargin{:});

addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = i_p.Results.I_folder;
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

gel_image_range = [Inf -Inf];
puncta_image_range = [Inf -Inf];
ranges = zeros(size(image_dirs), 1);
maxes = zeros(size(image_dirs), 1);
mins = zeros(size(image_dirs), 1);

for i = 1:size(image_dirs)
    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,'binary_shift.png')));
    
    gel_image = imread(fullfile(base_dir,image_dirs(i).name,'registered_gel.png'));
    if (min(gel_image(binary_shift)) < gel_image_range(1))
        gel_image_range(1) = min(gel_image(binary_shift));
    end
    if (max(gel_image(binary_shift)) > gel_image_range(2))
        gel_image_range(2) = max(gel_image(binary_shift));
    end
    
    ranges(i) = max(gel_image(binary_shift)) - min(gel_image(binary_shift));
    maxes(i) = max(gel_image(binary_shift));
    mins(i) = min(gel_image(binary_shift));
    
    puncta_image = imread(fullfile(base_dir,image_dirs(i).name,'registered_focal_image.png'));
    if (min(puncta_image(binary_shift)) < puncta_image_range(1))
        puncta_image_range(1) = min(puncta_image(binary_shift));
    end
    if (max(puncta_image(binary_shift)) > puncta_image_range(2))
        puncta_image_range(2) = max(puncta_image(binary_shift));
    end
    
end

for i = 1:size(image_dirs)
    csvwrite(fullfile(base_dir,image_dirs(i).name,'gel_image_range.csv'),gel_image_range)
    csvwrite(fullfile(base_dir,image_dirs(i).name,'puncta_image_range.csv'),puncta_image_range)
end

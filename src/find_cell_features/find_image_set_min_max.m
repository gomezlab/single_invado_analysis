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
filenames = add_filenames_to_struct(struct());
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

for i = 1:size(image_dirs,1)
    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,filenames.binary_shift)));
    
    gel_image = imread(fullfile(base_dir,image_dirs(i).name,filenames.gel));
    if (min(gel_image(binary_shift)) < gel_image_range(1))
        gel_image_range(1) = min(gel_image(binary_shift));
    end
    if (max(gel_image(binary_shift)) > gel_image_range(2))
        gel_image_range(2) = max(gel_image(binary_shift));
    end
    
    ranges(i) = max(gel_image(binary_shift)) - min(gel_image(binary_shift));
    maxes(i) = max(gel_image(binary_shift));
    mins(i) = min(gel_image(binary_shift));
    
    puncta_image = imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta));
    if (min(puncta_image(binary_shift)) < puncta_image_range(1))
        puncta_image_range(1) = min(puncta_image(binary_shift));
    end
    if (max(puncta_image(binary_shift)) > puncta_image_range(2))
        puncta_image_range(2) = max(puncta_image(binary_shift));
    end
    
    disp(['Done with ',num2str(i),'/',num2str(size(image_dirs,1))])
end

output_file = fullfile(base_dir,image_dirs(i).name,filenames.gel_range);
if(not(exist(fileparts(output_file),'dir')))
    mkdir(fileparts(output_file))
end
csvwrite(output_file,gel_image_range)

csvwrite(fullfile(base_dir,image_dirs(i).name,filenames.puncta_range),gel_image_range)
function find_gel_thresh(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'FIND_GEL_THRESH';

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

all_diff = [];

first_image = imread(fullfile(base_dir,image_dirs(1).name,'registered_gel.png'));
scale_factor = double(intmax(class(first_image)));
first_image  = double(first_image)/scale_factor;

for i = 1:size(image_dirs)
    gel_image = imread(fullfile(base_dir,image_dirs(i).name,'registered_gel.png'));
    scale_factor = double(intmax(class(gel_image)));
    gel_image  = double(gel_image)/scale_factor;
    
    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,'binary_shift.png')));
    
    diff_image = gel_image - first_image;
    diff_image = diff_image(binary_shift);
    assert(size(diff_image,1)*size(diff_image,2) == sum(sum(binary_shift)));
    
    all_diff = [all_diff diff_image(:)']; %#ok<AGROW>
end

threshold = mean(all_diff(:)) - 4*std(all_diff(:));

for i = 1:size(image_dirs)
    csvwrite(fullfile(base_dir,image_dirs(i).name,'gel_threshold.csv'),threshold)
end


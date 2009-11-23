function find_gel_thresh_hp(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_GEL_THRESH';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);

i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
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
    
    min_row = find(sum(binary_shift,2),1,'first');
    max_row = find(sum(binary_shift,2),1,'last');
    min_col = find(sum(binary_shift),1,'first');
    max_col = find(sum(binary_shift),1,'last');
    
    only_reg_first_image = first_image(min_row:max_row, min_col:max_col);
    only_reg_gel_image = gel_image(min_row:max_row, min_col:max_col);
    
    diff_image = only_reg_gel_image - only_reg_first_image;
    
    I_filt = fspecial('disk',i_p.Results.filter_size);
    
    blurred_image = imfilter(diff_image,I_filt,'same',mean(diff_image(:)));
    high_passed_image = diff_image - blurred_image;
    
    assert(size(high_passed_image,1)*size(high_passed_image,2) == sum(sum(binary_shift)));
    
    all_diff = [all_diff high_passed_image(:)']; %#ok<AGROW>
    
    if (i_p.Results.debug)
        disp(i)
    end
end

threshold = mean(all_diff(:)) - 3*std(all_diff(:));

for i = 1:size(image_dirs)
    csvwrite(fullfile(base_dir,image_dirs(i).name,'gel_threshold.csv'),threshold)
end


function find_puncta_thresh(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'find_puncta_thresh';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('image_file','registered_focal_image.png',@(x)exist(x,'file') == 2);
i_p.addParamValue('std_coeff',3,@(x)isnumeric(x))
i_p.addParamValue('output_file_name','puncta_threshold.csv',@(x)ischar(x))
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

all_filt = [];

for i = 1:size(image_dirs)
    image = imread(fullfile(base_dir,image_dirs(i).name,i_p.Results.image_file));
    scale_factor = double(intmax(class(image)));
    image  = double(image)/scale_factor;    

    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,'binary_shift.png')));
    min_row = find(sum(binary_shift,2),1,'first');
    max_row = find(sum(binary_shift,2),1,'last');
    min_col = find(sum(binary_shift),1,'first');
    max_col = find(sum(binary_shift),1,'last');
    
    only_reg_image = image(min_row:max_row, min_col:max_col);
    
    I_filt = fspecial('disk',i_p.Results.filter_size);
    blurred_image = imfilter(only_reg_image,I_filt,'same',mean(only_reg_image(:)));
    high_passed_image = only_reg_image - blurred_image;
    
    all_filt = [all_filt high_passed_image(:)']; %#ok<AGROW>
    if (i_p.Results.debug), disp(i); end
end

threshold = mean(all_filt(:)) + i_p.Results.std_coeff*std(all_filt(:));

for i = 1:size(image_dirs)
    csvwrite(fullfile(base_dir,image_dirs(i).name,i_p.Results.output_file_name),threshold)
end
function find_puncta_thresh(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'find_puncta_thresh';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('std_coeff',3,@(x)isnumeric(x))
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_folder,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(i_p.Results.I_folder,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

all_filt = [];

for i = 1:size(image_dirs)
    image = imread(fullfile(base_dir,image_dirs(i).name,filenames.gel));

    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,filenames.binary_shift)));
    
    only_reg_image = remove_region_outside_registered(image,binary_shift);

    I_filt = fspecial('disk',11);
    blurred_image = imfilter(only_reg_image,I_filt,'same',mean(only_reg_image(:)));
    high_passed_image = only_reg_image - blurred_image;
    
    all_filt = [all_filt high_passed_image(:)']; %#ok<AGROW>
    disp(['Done with ',num2str(i),'/',num2str(size(image_dirs,1))]);
end
all_filt = double(all_filt);

threshold = mean(all_filt(:)) + i_p.Results.std_coeff*std(all_filt(:));

csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold),threshold)

toc;
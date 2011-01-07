function flat_field_correct_puncta(base_dir, varargin)

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

flat_field = double(imread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_median)));

for i = 1:size(image_dirs)
    if (exist(fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.puncta_filename]),'file'))
        continue;
    end
    
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta_filename)));
    
    puncta_image = puncta_image - flat_field;
    puncta_image(puncta_image < 0) = 0;
    
    movefile(fullfile(base_dir,image_dirs(i).name,filenames.puncta_filename), ... 
        fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.puncta_filename]));
    
    imwrite(uint16(puncta_image),fullfile(base_dir,image_dirs(i).name,filenames.puncta_filename),'bitdepth',16)
end
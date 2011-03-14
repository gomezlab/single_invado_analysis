function cascade_registration_transformations(I_folder, varargin)
% CASCADE_REGISTRATION_TRANSFORMATIONS    takes a set of adjacent image
%       registrations and maintains a running total of the x and y
%       transformations, writing these new cascaded transformations to the
%       same folder  
%
%   cascade_registration_transformations(I_dir) reads in all the
%   transformation matrices in image folders (I_dir) starting at image 2
%   and keeps a running total of the x and y transformations, writing those
%   cascaded transformations back to the same transformation matrices
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'CASCADE_REGISTRATION_TRANSFORMATIONS';

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

x_transform = 0;
y_transform = 0;

for i = 1:size(image_dirs)
    transform_file = fullfile(base_dir,image_dirs(i).name,filenames.affine_matrix);
    
    assert(exist(transform_file,'file') == 2)
    transform_mat = csvread(transform_file);
    assert(all(all(transform_mat(1:2,1:2) == [1 0; 0 1])))
    
    x_transform = x_transform + transform_mat(3,1);
    y_transform = y_transform + transform_mat(3,2);
    if (i_p.Results.debug)
        disp(sprintf('%d,%d',x_transform, y_transform)) %#ok<DSPS>
    end
    new_transform_mat = transform_mat;
    new_transform_mat(3,1) = x_transform;
    new_transform_mat(3,2) = y_transform;
    
    new_transform_file = fullfile(base_dir,image_dirs(i).name,filenames.affine_matrix_cascade);
    
    csvwrite(new_transform_file, new_transform_mat);
end

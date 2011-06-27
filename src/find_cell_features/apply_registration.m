function apply_registration(I_folder,varargin)
% APPLY_REGISTRATION    Applies the given transformation matrix to an
%                       image or if no matrix is provided, applies no
%                       tranformation
%
%   apply_registration(I1,OPTIONS) Applies the default registration (no
%   change) to the image (I1)
%
%   Options:
%
%       -transform_file: csv file containing the affine tranformation to be
%        applied to the image
%       -output_dir: folder used to hold all the results, defaults to the
%        same folder as the image file, 'I'
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'APPLY_REGISTRATION';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_folder,varargin{:});

%Add the folder with all the scripts used in this master program
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

for i = 1:size(image_dirs,1)
    puncta_image  = imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta));
    gel_image  = imread(fullfile(base_dir,image_dirs(i).name,filenames.gel));
    transform_matrix = csvread(fullfile(base_dir,image_dirs(i).name,filenames.affine_matrix_cascade));
    transform = maketform('affine', transform_matrix);
    
    binary_image = ones(size(puncta_image));
    
    binary_shift = imtransform(binary_image, transform, 'XData',[1 size(puncta_image,2)], 'YData', [1 size(puncta_image,1)]);
    puncta_shift = imtransform(puncta_image, transform, 'XData',[1 size(puncta_image,2)], 'YData', [1 size(puncta_image,1)]);
    gel_shift = imtransform(gel_image, transform, 'XData',[1 size(puncta_image,2)], 'YData', [1 size(puncta_image,1)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imwrite(binary_shift,fullfile(base_dir,image_dirs(i).name,filenames.binary_shift));
    imwrite(gel_shift,fullfile(base_dir,image_dirs(i).name,filenames.gel));
    imwrite(puncta_shift,fullfile(base_dir,image_dirs(i).name,filenames.puncta));
    
    disp(['Done with ',num2str(i),'/',num2str(size(image_dirs,1))])
end

toc;
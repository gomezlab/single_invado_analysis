function apply_registration(I_file,varargin)
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

i_p = inputParser;
i_p.FunctionName = 'APPLY_REGISTRATION';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.addParamValue('transform_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%read in and normalize the input focal adhesion image
orig_image  = imread(I_file);
scale_factor = double(intmax(class(orig_image)));
orig_image  = double(orig_image)/scale_factor;

if (isempty(strmatch('transform_file',i_p. UsingDefaults)))
    transform_matrix = csvread(i_p.Results.transform_file);
else
    transform_matrix = [1 0; 0 1; 0 0];
end

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transform = maketform('affine', transform_matrix);

binary_image = ones(size(orig_image));
binary_shift = imtransform(binary_image, transform, 'XData',[1 size(orig_image,2)], 'YData', [1 size(orig_image,1)]);
image_shift = imtransform(orig_image, transform, 'XData',[1 size(orig_image,2)], 'YData', [1 size(orig_image,1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[path, name, ext, version] = fileparts(i_p.Results.I_file);

imwrite(binary_shift,fullfile(i_p.Results.output_dir, 'binary_shift.png'));
imwrite(image_shift,fullfile(i_p.Results.output_dir, ['registered_' name '.png']),'bitdepth',16);
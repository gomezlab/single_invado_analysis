function [varargout] = register_gel_images(I_file,first_image,varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,OPTIONS) locate the focal adhesions in image
%   file, 'I'
%
%   Options:
%
%       -cell_mask: file which contains the cell mask
%       -filter_size: size of the averaging filter to use, defaults to 23
%       -filter_thresh: threshold used to identify focal adhesions in the
%        average filtered image, defaults to 0.1
%       -output_dir: folder used to hold all the results, defaults to the
%        same folder as the image file, 'I'
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'REGISTER_GEL_IMAGES';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('first_image',@(x)exist(x,'file') == 2);
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,first_image,varargin{:});

%read in and normalize the input focal adhesion image
gel_image  = imread(I_file);
scale_factor = double(intmax(class(gel_image)));
gel_image  = double(gel_image)/scale_factor;

%read in and normalize the first image
first_image  = imread(first_image);
scale_factor = double(intmax(class(first_image)));
first_image  = double(first_image)/scale_factor;

    
%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row_shifts,col_shifts] = meshgrid(-100:10:100);
ms_diff = zeros(size(row_shifts));

for i = 1:size(row_shifts,1)
    disp(i)
    for j = 1:size(col_shifts,2)
        this_row_shift = row_shifts(i,j);
        this_col_shift = col_shifts(i,j);

        transform = maketform('affine',[cos(0) sin(0);-sin(0) cos(0); this_row_shift this_col_shift]);
        
        binary_image = ones(size(gel_image));
        binary_shift = imtransform(binary_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
        gel_shift = imtransform(gel_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
        
        ms_diff(i,j) = sum(sum((first_image.*binary_shift - gel_shift.*binary_shift).^2))/sum(sum(binary_shift));
    end
end

best_index = find(ms_diff == min(min(ms_diff)),1)
transform_matrix = [cos(0) sin(0);-sin(0) cos(0); row_shifts(best_index) col_shifts(best_index)];
transform = maketform('affine', transform_matrix);

binary_image = ones(size(gel_image));
binary_shift = imtransform(binary_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
gel_shift = imtransform(gel_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(binary_shift,fullfile(i_p.Results.output_dir, 'binary_shift.png'));
imwrite(gel_shift,fullfile(i_p.Results.output_dir, 'registered_gel.png'));
csvwrite(fullfile(i_p.Results.output_dir, 'affine_matrix.csv'), transform_matrix)
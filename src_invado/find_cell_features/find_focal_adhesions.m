function find_focal_adhesions(I_file,varargin)
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
i_p.FunctionName = 'FIND_FOCAL_ADHESIONS';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

i_p.addParamValue('min_size',0.56,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('binary_shift_file',fullfile(fileparts(I_file),'binary_shift.png'),@(x)exist(x,'file')==2);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',0.2,@isnumeric);
i_p.addParamValue('scale_filter_thresh',0,@(x)islogical(x) || (isnumeric(x) && (x == 1 || x == 0)));
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%read in binary shift file
binary_shift = logical(imread(i_p.Results.binary_shift_file));

%read in and normalize the input focal adhesion image
focal_image  = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

if (i_p.Results.scale_filter_thresh) 
    filter_thresh = i_p.Results.filter_thresh * (max(focal_image(:)) - min(focal_image(:)));
else
    filter_thresh = i_p.Results.filter_thresh;
end
    
%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(focal_image(binary_shift),I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image(binary_shift) - blurred_image;
threshed_image = logical(im2bw(high_passed_image,filter_thresh));

full_size_thresh = zeros(size(focal_image));
full_size_thresh(binary_shift) = threshed_image;
threshed_image = full_size_thresh;

puncta = bwlabel(threshed_image,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(double(puncta)/2^16,fullfile(i_p.Results.output_dir, 'puncta_labeled.png'),'bitdepth',16);
imwrite(threshed_image,fullfile(i_p.Results.output_dir, 'puncta_binary.png'));

addpath(genpath('..'))

only_registered = focal_image(binary_shift);

scaled_image = focal_image;
scaled_image = scaled_image - min(only_registered(:));
scaled_image = scaled_image .* (1/max(only_registered(:)));
scaled_image(not(binary_shift)) = 0;

imwrite(create_highlighted_image(scaled_image, threshed_image),fullfile(i_p.Results.output_dir, 'puncta_highlight.png')); 
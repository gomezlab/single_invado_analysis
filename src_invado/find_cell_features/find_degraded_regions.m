function find_degraded_regions(I_file,first_image,varargin)
% FIND_DEGRADED_REGIONS    locates the regions of degradation (decrease in
%                          fluorescence intensity) between two images
%                         
%   find_degraded_regions(I1,I2,OPTIONS) searches for regions in the first
%   image (I1) where the fluorescence has decreased as compared to the
%   second image (I2)
%
%   Options:
%
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
i_p.FunctionName = 'FIND_DEGRADED_REGIONS';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

i_p.addRequired('first_image',@(x)exist(x,'file') == 2);
i_p.addParamValue('registration_binary',@(x)exist(x,'file') == 2);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',-0.075,@isnumeric);
i_p.addParamValue('scale_filter_thresh',0,@(x)islogical(x) || (isnumeric(x) && (x == 1 || x == 0)));
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

%read in and normalize the registration binary image
if (isempty(strmatch('registration_binary',i_p.UsingDefaults)))
    binary_shift = logical(imread(i_p.Results.registration_binary));
    first_image = first_image .* binary_shift;
end

if (i_p.Results.scale_filter_thresh)
    filter_thresh = i_p.Results.filter_thresh * (max(gel_image(:)) - min(gel_image(:)));
else
    filter_thresh = i_p.Results.filter_thresh;
end
    
%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
% blurred_image = imfilter(gel_image,I_filt,'same',mean(gel_image(:)));
% high_passed_image = gel_image - blurred_image;
% threshed_gel = high_passed_image < filter_thresh;
% 
% I_filt = fspecial('disk',i_p.Results.filter_size);
% blurred_image = imfilter(first_image,I_filt,'same',mean(first_image(:)));
% high_passed_image = first_image - blurred_image;
% threshed_first = high_passed_image < filter_thresh;

diff_image = gel_image - first_image;
threshed_diff = diff_image < i_p.Results.filter_thresh;



blurred_image = imfilter(diff_image,I_filt,'same',mean(diff_image(:)));
high_passed_image = diff_image - blurred_image;
threshed_diff_filt = high_passed_image < i_p.Results.filter_thresh;
threshed_diff_filt = logical(threshed_diff_filt .* binary_shift);

addpath(genpath('..'))

both_diff = zeros(size(gel_image));
both_diff(threshed_diff) = 1;
both_diff(threshed_diff_filt) = 2;

if (i_p.Results.debug)
    imshow(create_highlighted_image(gel_image,both_diff, 'color_map', [1 0 0; 0 1 0]));
    figure;
    imshow(create_highlighted_image(gel_image,threshed_diff));
    figure;
    imshow(create_highlighted_image(gel_image,threshed_diff_filt));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imwrite(threshed_diff,fullfile(i_p.Results.output_dir, 'degradation_binary.png')); 

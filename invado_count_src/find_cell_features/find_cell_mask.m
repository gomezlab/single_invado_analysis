function [varargout] = find_cell_mask(I_file,out_file,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('out_file',@(x)isempty(fileparts(x)) == 1 || exist(fileparts(x),'dir') == 7);

i_p.addParamValue('binary_shift_file','',@(x)exist(x,'file') == 2);

i_p.parse(I_file,out_file,varargin{:});

mask_image = imread(I_file);
scale_factor = double(intmax(class(mask_image)));
mask_image   = double(mask_image)/scale_factor;

pixel_values = [];
if (not(isempty(strmatch('binary_shift_file',i_p.UsingDefaults))))
    pixel_values = mask_image(:);
else
    binary_shift = logical(imread(i_p.Results.binary_shift_file));
    pixel_values = mask_image(binary_shift);
    assert(length(pixel_values) == sum(sum(binary_shift)));
end


puncta_min_max = csvread(fullfile(fileparts(I_file),'puncta_image_range.csv'));
puncta_min_max = puncta_min_max/scale_factor;

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Threshold identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorted_mask_pixels = sort(pixel_values);
% sorted_mask_pixels(1:0.05*round(length(sorted_mask_pixels))) = 0;

%when there are very few unique pixel values, having a large number of bins
%causes the extrema values found to be fairly random, fixing with a simple
%check for the number of uniques
if (length(unique(sorted_mask_pixels)) < 1000)
    [heights, intensity] = hist(sorted_mask_pixels,100);
else
    [heights, intensity] = hist(sorted_mask_pixels,1000);
end

if (length(unique(sorted_mask_pixels)) < 1000)
    smoothed_heights = smooth(heights,0.15,'loess');
else
    smoothed_heights = smooth(heights,0.05,'loess');
end

[zmax,imax,zmin,imin]= extrema(smoothed_heights);

%diagnostic plot for when things don't go quite right
% plot(intensity, smoothed_heights)
% hold on;
% plot(intensity(imax), zmax,'gx')
% plot(intensity(imin), zmin,'ro')

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin
sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));

threshed_mask = im2bw(mask_image, intensity(imin(min_index)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Mask Cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshed_mask = imfill(threshed_mask,'holes');

connected_areas = bwlabel(threshed_mask,8);
region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>

%filter out connected regions smaller than 100 pixels
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 100));
connected_areas = bwlabel(threshed_mask,8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output Image Creation/Writing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cell mask highlighting
normalized_image = mask_image - puncta_min_max(1);
normalized_image = normalized_image / (puncta_min_max(2) - puncta_min_max(1));

imwrite(create_highlighted_image(normalized_image,bwperim(threshed_mask)), ...
    fullfile(fileparts(out_file),'highlighted_mask.png'))

%binary cell mask
imwrite(threshed_mask, out_file)

%labeled cell mask
imwrite(double(connected_areas)/2^16, ...
    fullfile(fileparts(out_file),'cell_mask_labeled.png'), ...
    'bitdepth',16)

if (nargout >= 1)
    varargout{1} = threshed_mask;
end

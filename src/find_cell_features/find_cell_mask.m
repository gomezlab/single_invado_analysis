function [varargout] = find_cell_mask(I_file,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in images and properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_folder = fileparts(I_file);

mask_image = double(imread(I_file));

binary_shift = logical(imread(fullfile(I_folder,filenames.binary_shift)));
only_reg_pixels = mask_image(binary_shift);
assert(length(only_reg_pixels) == sum(sum(binary_shift)));

puncta_min_max = csvread(fullfile(I_folder,filenames.puncta_range));


%%Threshold identification
sorted_mask_pixels = sort(only_reg_pixels);
% sorted_mask_pixels(1:0.05*round(length(sorted_mask_pixels))) = 0;

%when there are very few unique pixel values, having a large number of bins
%causes the extrema values found to be fairly random, fixing with a simple
%check for the number of uniques
% if (length(unique(sorted_mask_pixels)) < 1000)
%     [heights, intensity] = hist(sorted_mask_pixels,100);
% else
%     [heights, intensity] = hist(sorted_mask_pixels,1000);
% end

[heights, intensity] = hist(sorted_mask_pixels,length(unique(sorted_mask_pixels))/5);

% similar problem with few pixel values
% if (length(unique(sorted_mask_pixels)) < 1000)
%     smoothed_heights = smooth(heights,0.15,'loess');
% else
%     smoothed_heights = smooth(heights,0.05,'loess');
% end

smoothed_heights = smooth(heights,0.10,'loess');

[zmax,imax,zmin,imin]= extrema(smoothed_heights);

if(i_p.Results.debug)
    plot(intensity, smoothed_heights)
    hold on;
    plot(intensity(imax), zmax,'gx')
    plot(intensity(imin), zmin,'ro')
end

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin

sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));

threshed_mask = mask_image > intensity(imin(min_index));

%%Mask Cleanup
threshed_mask = imfill(threshed_mask,'holes');

connected_areas = bwlabel(threshed_mask);
region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>

%filter out connected regions smaller than 10000 pixels
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 10000));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image writing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imwrite(threshed_mask, fullfile(I_folder,filenames.cell_mask));

%scale the input image and write out a highlighted image showing the
%detected cell edge
normalized_image = (mask_image - puncta_min_max(1))/range(puncta_min_max);

imwrite(create_highlighted_image(normalized_image,bwperim(threshed_mask)),fullfile(I_folder,filenames.cell_mask_highlight))

if (nargout >= 1)
    varargout{1} = threshed_mask;
end

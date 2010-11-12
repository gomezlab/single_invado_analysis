function [varargout] = find_global_cell_mask(base_dir,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);

i_p.parse(base_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

image_folders = dir(base_dir);
image_folders = image_folders(3:end);

mask_images = cell(0);
binary_shifts = cell(0);

pixel_values = []; %#ok<NASGU>
for i=1:length(image_folders)
    mask_images{i} = imread(fullfile(base_dir,image_folders(i).name,filenames.puncta_filename));
    puncta_min_max = csvread(fullfile(base_dir,image_folders(i).name,filenames.puncta_range_file));
    mask_images{i} = (double(mask_images{i}) - puncta_min_max(1)) / range(puncta_min_max);
    
    binary_shifts{i} = logical(imread(fullfile(base_dir,image_folders(i).name,filenames.binary_shift_filename)));
    mask_images{i}(not(binary_shifts{i})) = 0;
end

pixel_values = [];
for i=1:length(mask_images)
    pixel_values_temp = mask_images{i}(binary_shifts{i});
    pixel_values = cat(1,pixel_values,pixel_values_temp(:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Threshold identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorted_mask_pixels = double(sort(pixel_values));
% sorted_mask_pixels(1:0.05*round(length(sorted_mask_pixels))) = 0;
sorted_mask_pixels = sorted_mask_pixels(0.05*round(length(sorted_mask_pixels)):end);

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
plot(intensity, smoothed_heights)
hold on;
plot(intensity(imax), zmax,'gx')
plot(intensity(imin), zmin,'ro')

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin
sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));
cell_mask_threshold = intensity(imin(min_index));

cell_masks = cell(0);
labeled_areas = cell(0);
for i=1:length(mask_images)
    cell_masks{i} = im2bw(mask_images{i}, cell_mask_threshold);
    cell_masks{i} = clean_mask_image(cell_masks{i});
    
    if (sum(sum(cell_masks{i})) > 0)
        1;
    end
    
end

for i=1:length(image_folders)
    imwrite(cell_masks{i},fullfile(base_dir,image_folders(i).name,filenames.cell_mask_filename));

    connected_areas = bwlabel(cell_masks{i},8);
    
    imwrite(create_highlighted_image(mask_images{i},bwperim(cell_masks{i})), ...
        fullfile(base_dir,image_folders(i).name,filenames.highlighted_cell_mask_filename))

    imwrite(double(connected_areas)/2^16, ...
        fullfile(base_dir,image_folders(i).name,filenames.labeled_cell_mask_filename), ...
        'bitdepth',16)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output Image Creation/Writing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cell mask highlighting and labeled cell mask output, if we are working
%with the registered image
% if (isempty(strmatch('binary_shift_file',i_p.UsingDefaults)))
%     normalized_image = mask_image - puncta_min_max(1);
%     normalized_image = normalized_image / (puncta_min_max(2) - puncta_min_max(1));
%     
%     imwrite(create_highlighted_image(normalized_image,bwperim(threshed_mask)), ...
%         fullfile(fileparts(out_file),'highlighted_mask.png'))
%     
%     imwrite(double(connected_areas)/2^16, ...
%         fullfile(fileparts(out_file),'cell_mask_labeled.png'), ...
%         'bitdepth',16)
% end
% 
% %binary cell mask
% imwrite(threshed_mask, out_file)
% 
% if (nargout >= 1)
%     varargout{1} = threshed_mask;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [threshed_mask] = clean_mask_image(threshed_mask)

threshed_mask = imfill(threshed_mask,'holes');

connected_areas = bwlabel(threshed_mask,8);
region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>

%filter out connected regions smaller than 100 pixels
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 100));

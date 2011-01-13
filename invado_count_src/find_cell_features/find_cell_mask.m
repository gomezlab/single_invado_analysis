function find_cell_mask(exp_dir,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('min_cell_area',1500,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

for i_num = 1:size(image_dirs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reading in current image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta_filename)));
    puncta_min_max = csvread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta_range_file));
    pixel_values = puncta_image(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reading in prior connected areas
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (i_num == 1)
        prior_connected_areas = [];
    else
        prior_connected_areas = double(imread(fullfile(base_dir,image_dirs(i_num-1).name,filenames.labeled_cell_mask_filename)));
        1;
    end
    
    thresh = find_mask_threshold(pixel_values,i_p);    
    threshed_mask = puncta_image > thresh;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Mask Cleanup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    connected_areas = filter_mask(threshed_mask, puncta_image, i_p);
    threshed_mask = connected_areas > 0;
    
    if (not(isempty(prior_connected_areas)))
        connected_areas = filter_on_overlap(puncta_image,connected_areas,prior_connected_areas);
        threshed_mask = connected_areas > 0;
        connected_areas = filter_mask(threshed_mask, puncta_image, i_p);
        threshed_mask = connected_areas > 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Output Image Creation/Writing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %cell mask highlighting and labeled cell mask output, if we are working
    %with the registered image
    normalized_image = puncta_image - puncta_min_max(1);
    normalized_image = normalized_image / (puncta_min_max(2) - puncta_min_max(1));
    
    imwrite(create_highlighted_image(normalized_image,bwperim(threshed_mask),'color_map',[0,1,0]), ...
        fullfile(base_dir,image_dirs(i_num).name,filenames.highlighted_cell_mask_filename))
    
    imwrite(double(connected_areas)/2^16, ...
        fullfile(base_dir,image_dirs(i_num).name,filenames.labeled_cell_mask_filename), ...
        'bitdepth',16)

    %binary cell mask
    
    imwrite(threshed_mask, fullfile(base_dir,image_dirs(i_num).name,filenames.cell_mask_filename))
    if (mod(i_num,10)==0)
        disp(['Done processing image number: ',num2str(i_num)])
    end
    disp(['Done processing image number: ',num2str(i_num)])
end
toc;

function thresh = find_mask_threshold(pixel_values,i_p)

sorted_mask_pixels = sort(pixel_values);
sorted_mask_pixels = sorted_mask_pixels(round(0.05*length(sorted_mask_pixels)):end);
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

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin
sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));

%diagnostic plot for when things don't go quite right
if (i_p.Results.debug)
    plot(intensity, smoothed_heights)
    hold on;
    plot(intensity(imax), zmax,'gx')
    plot(intensity(imin), zmin,'ro')
    plot(intensity(imin(min_index)), zmin(min_index),'k*','MarkerSize',16);
end

thresh = intensity(imin(min_index));

function connected_areas = filter_mask(threshed_mask, puncta_image, i_p)

threshed_mask = imfill(threshed_mask,'holes');

connected_areas = bwlabel(threshed_mask,8);
region_props = regionprops(connected_areas,puncta_image, 'Area','MeanIntensity'); %#ok<MRPBW>

%filter out connected regions smaller than the min cell area
% size_filter = [region_props.Area] > i_p.Results.min_cell_area &  ...
%     [region_props.Area] < 10000;

size_filter = [region_props.Area] > i_p.Results.min_cell_area;

intensity_filter = [region_props.MeanIntensity] > 0.002;

threshed_mask = ismember(connected_areas, find(size_filter & intensity_filter));

connected_areas = bwlabel(threshed_mask,8);

function final_connected_areas = filter_on_overlap(puncta_image,connected_areas,prior_connected_areas)

final_connected_areas = zeros(size(puncta_image));

for (i=1:max(connected_areas(:)))
    this_connected_area = connected_areas == i;
    con_props = regionprops(this_connected_area,'Area');
    
    seeds = prior_connected_areas > 0 & this_connected_area;
    seeds_labeled = bwlabel(seeds,8);
    seeds_props = regionprops(seeds_labeled,'Area');
%     if (max(seeds_labeled(:)) <= 1 || (sum([seeds_props.Area]) > 0.9*sum([con_props.Area])))
     if (max(seeds_labeled(:)) <= 1)
        final_connected_areas(this_connected_area) = max(final_connected_areas(:)) + 1;
    else
        puncta_image_invert = puncta_image*-1+max(puncta_image(:));
        
        no_cells = not(this_connected_area);
        
        minned_image = imimposemin(puncta_image_invert,seeds,8);
        
        water_out = watershed(minned_image);
        water_out(no_cells) = 0;
        
        for (i=1:max(water_out(:)))
            final_connected_areas(water_out == i) = max(final_connected_areas(:)) + 1;
        end    
    end
    
end

cell_nums = unique(final_connected_areas);
assert(cell_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(cell_nums)
    final_connected_areas(final_connected_areas == cell_nums(i)) = i - 1;
end



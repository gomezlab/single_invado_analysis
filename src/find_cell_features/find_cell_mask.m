function find_cell_mask(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholds = NaN*ones(size(image_dirs,1),1);

for i_num = 1:size(image_dirs)
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta)));

    thresholds(i_num) = find_image_threshold(puncta_image(:),i_p.Results.debug);
end

overall_threshold = nanmedian(thresholds(:));
puncta_min_max = csvread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_range));

for i_num = 1:size(image_dirs)
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta)));
    
    threshed_mask = puncta_image > overall_threshold;
    
    %%Mask Cleanup
    threshed_mask = imfill(threshed_mask,'holes');
    
    connected_areas = bwlabel(threshed_mask);
    region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>
    
    %filter out connected regions smaller than 10000 pixels
    threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 10000));
    
    imwrite(threshed_mask, fullfile(base_dir,image_dirs(i_num).name,filenames.cell_mask));
    if (mod(i_num,10)==0)
        disp(['Done processing image number: ',num2str(i_num)])
    end
    
    if (i_p.Results.debug)
        puncta_norm = (puncta_image - puncta_min_max(1))/range(puncta_min_max);
        edge_highlight = create_highlighted_image(puncta_norm,bwperim(threshed_mask));
        imwrite(edge_highlight, fullfile(base_dir,image_dirs(i_num).name,'edge_highlight.png'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold = find_image_threshold(puncta_image_pixels,debug)

sorted_mask_pixels = sort(puncta_image_pixels);
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

smoothed_heights = smooth(heights,0.10,'loess');

[zmax,imax,zmin,imin]= extrema(smoothed_heights);

if(debug)
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
try
    min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
    assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));
    threshold = intensity(imin(min_index));
catch %#ok<CTCH>
    threshold = NaN;
end

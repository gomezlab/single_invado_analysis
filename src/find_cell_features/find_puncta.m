function find_puncta(I_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.addParamValue('min_puncta_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'))
filenames = add_filenames_to_struct(struct());
I_folder = fileparts(I_file);

%read in binary shift file
binary_shift = logical(imread(fullfile(I_folder,filenames.binary_shift)));

cell_mask = logical(imread(fullfile(I_folder,filenames.cell_mask)));

%read in and remove regions outside registered region
focal_image  = double(imread(I_file));
only_reg_focal_image = remove_region_outside_registered(focal_image,binary_shift);

%read in the global min max file
min_max = csvread(fullfile(I_folder,filenames.puncta_range));

filter_thresh = csvread(fullfile(I_folder,filenames.puncta_threshold));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(only_reg_focal_image,I_filt,'same',mean(only_reg_focal_image(:)));
high_passed_image = only_reg_focal_image - blurred_image;

threshed_image = find_threshed_image(high_passed_image,filter_thresh);

%identify and remove objects on the immediate edge of the image
threshed_image = remove_edge_objects(threshed_image);

%place the thresholded image back in place
threshed_temp = zeros(size(focal_image));
threshed_temp(binary_shift) = threshed_image;
threshed_image = threshed_temp;

if (i_p.Results.min_puncta_size > 1)
    labeled_thresh = bwlabel(threshed_image,4);
    
    props = regionprops(labeled_thresh,'Area'); %#ok<MRPBW>
    labeled_thresh = ismember(labeled_thresh, find([props.Area] >= i_p.Results.min_puncta_size));
    
    threshed_image = labeled_thresh > 0;
end

puncta = bwlabel(threshed_image,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove objects outside mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:max(puncta(:))
    assert(any(any(puncta == i)), 'Error: can''t find ad number %d', i);
    this_puncta = zeros(size(puncta));
    this_puncta(puncta == i) = 1;
    if (sum(sum(this_puncta & cell_mask)) == 0)
        puncta(puncta == i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and fill holes in single puncta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
props = regionprops(puncta,'Area');
large_ad_nums = find([props.Area] >= 4);
for this_num = large_ad_nums
    %first make a binary image of the current adhesion and then run imfill
    %to fill any holes present    
    this_puncta = puncta == this_num;
    filled_ad = imfill(this_puncta,'holes');
    
    puncta(logical(filled_ad)) = this_num;
    if (i_p.Results.debug && mod(this_num,50)==0)
        disp(['Done filling holes in ',num2str(this_num), '/', num2str(length(ad_nums))]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Placeholder image for empties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (sum(sum(puncta)) == 0)
    puncta(end,end) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumber objects to be sequential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(puncta);
assert(ad_nums(1) == 0, 'Background pixels not found after building puncta label matrix')
for i = 2:length(ad_nums)
    puncta(puncta == ad_nums(i)) = i - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build puncta perimeters image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
puncta_perim = zeros(size(puncta));
for i = 1:max(puncta(:))
    assert(any(any(puncta == i)), 'Error: can''t find ad number %d', i);
    this_puncta = zeros(size(puncta));
    this_puncta(puncta == i) = 1;
    puncta_perim(bwperim(this_puncta)) = i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(double(puncta)/2^16,fullfile(I_folder, filenames.objects),'bitdepth',16);
imwrite(double(puncta_perim)/2^16,fullfile(I_folder, filenames.objects_perim),'bitdepth',16);
imwrite(im2bw(puncta),fullfile(I_folder, filenames.objects_binary));

scaled_image = (focal_image - min_max(1))/(range(min_max));
scaled_image(not(binary_shift)) = 0;
highlighted_image = create_highlighted_image(scaled_image, im2bw(puncta),'color_map',[1 0 0]);
highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask),'color_map',[0 1 0]);

imwrite(highlighted_image, fullfile(I_folder, filenames.objects_highlight));
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshed_image = find_threshed_image(high_passed_image, filter_thresh)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_objects(high_threshed_image);
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    if (overlap_labels(1) == 0)
        overlap_labels = overlap_labels(2:end);
    end
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
end
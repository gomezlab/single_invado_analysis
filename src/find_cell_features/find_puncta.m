function find_puncta(I_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

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

%read in and normalize the input focal adhesion image
focal_image  = double(imread(I_file));
only_reg_focal_image = remove_region_outside_registered(focal_image,binary_shift);

%read in the global min max file, if present, will be used in the last part
%to output a normalized copy of the image with the puncta highlighted
min_max = csvread(fullfile(I_folder,filenames.puncta_range));

filter_thresh = csvread(fullfile(I_folder,filenames.puncta_threshold));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(only_reg_focal_image,I_filt,'same',mean(only_reg_focal_image(:)));
high_passed_image = only_reg_focal_image - blurred_image;
threshed_image = high_passed_image > filter_thresh;

%identify and remove objects on the immediate edge of the image
threshed_image = remove_edge_objects(threshed_image);

%place the thresholded image back in place
threshed_temp = zeros(size(focal_image));
threshed_temp(binary_shift) = threshed_image;
threshed_image = threshed_temp;

puncta = bwlabel(threshed_image,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove objects outside mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:max(puncta(:))
    assert(any(any(puncta == i)), 'Error: can''t find ad number %d', i);
    this_ad = zeros(size(puncta));
    this_ad(puncta == i) = 1;
    if (sum(sum(this_ad & cell_mask)) == 0)
        puncta(puncta == i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumber objects to be sequential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(puncta);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    puncta(puncta == ad_nums(i)) = i - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build adhesion perimeters image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
puncta_perim = zeros(size(puncta));
for i = 1:max(puncta(:))
    assert(any(any(puncta == i)), 'Error: can''t find ad number %d', i);
    this_ad = zeros(size(puncta));
    this_ad(puncta == i) = 1;
    puncta_perim(bwperim(this_ad)) = i;
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
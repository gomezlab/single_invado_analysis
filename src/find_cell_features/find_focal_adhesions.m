function [varargout] = find_focal_adhesions(I_file,varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,OPTIONS) locate the focal adhesions in image
%   file, 'I', and collects the properties of the focal adhesions
%
%   Options:
%
%       -cell_mask: file which contains the cell mask, defaults to not
%        present
%       -filter_size: size of the MATLAB average filter to use, defaults to
%        23
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

i_p.addOptional('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addOptional('min_size',0.56,@(x)isnumeric(x) && x > 1);
i_p.addOptional('pixel_size',0.215051,@isnumeric);
i_p.addOptional('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addOptional('filter_thresh',0.1,@isnumeric);
i_p.addOptional('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addOptional('output_file', 'adhesions.png', @ischar);
i_p.addOptional('output_file_perim', 'adhesions_perim.png', @ischar);
i_p.addOptional('output_file_binary', 'adhesions_binary.png', @ischar);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image  = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;

ad_zamir = find_ad_zamir(high_passed_image,i_p);
if (exist('cell_mask','var'))
    ad_zamir = find_in_cell_ads(ad_zamir,cell_mask);
    ad_nums = unique(ad_zamir);
    assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix');
    for i = 2:length(ad_nums)
        ad_zamir(ad_zamir == ad_nums(i)) = i - 1;
    end
end

ad_zamir_perim = zeros(size(ad_zamir));
for i = 1:max(ad_zamir(:))
    assert(any(any(ad_zamir == i)), 'Error: can''t find ad number %d', i);

    this_ad = zeros(size(ad_zamir));
    this_ad(ad_zamir == i) = 1;
    ad_zamir_perim(bwperim(this_ad)) = i;
end

if (exist('cell_mask','var'))
    adhesion_properties = collect_adhesion_properties(ad_zamir,focal_image,'cell_mask',cell_mask);
else
    adhesion_properties = collect_adhesion_properties(ad_zamir,focal_image);
end

%write the results to files
imwrite(double(ad_zamir)/2^16,fullfile(i_p.Results.output_dir, i_p.Results.output_file),'bitdepth',16);
imwrite(double(ad_zamir_perim)/2^16,fullfile(i_p.Results.output_dir, i_p.Results.output_file_perim),'bitdepth',16);
imwrite(im2bw(ad_zamir,0),fullfile(i_p.Results.output_dir, i_p.Results.output_file_binary));
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));

if (nargout > 0)
    varargout{1} = struct('adhesions',adhesions,'ad_zamir',ad_zamir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered_adhesions = find_in_cell_ads(ads,cm)
% COLLECT_ADHESION_PROPERTIES    using the identified adhesions, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                adhesions
%
%   ad_p = collect_adhesion_properties(ad_I,c_m,orig_I) collects the
%   properties of the adhesions identified in the binary image 'ad_I',
%   using the cell mask in 'c_m' and the original focal image data in
%   'orig_I', returning a structure 'ad_p' containing properties
%
%   Properties Collected:
%       -all of the properties collected by regioprops(...,'all')
%       -the distance of each adhesion's centroid from the nearest cell
%        edge
%       -the average and variance of the normalized fluorescence signal
%        within each adhesion

filtered_adhesions = zeros(size(ads,1),size(ads,2));

for i = 1:max(ads(:))
    this_adhesion = false(size(ads,1),size(ads,2));
    this_adhesion(ads == i) = true;
    overlap = and(this_adhesion,cm);

    if (sum(overlap(:)) >= 1)
        filtered_adhesions(this_adhesion) = i;
    end
end

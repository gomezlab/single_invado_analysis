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
%       -filt_size: size of the MATLAB average filter to use, defaults to
%        23
%       -min_intensity: threshold used to identify focal adhesions in the
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

i_p.addParamValue('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addParamValue('filt_size',23,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('min_intensity',0.1,@isnumeric);
i_p.addParamValue('output_dir',fileparts(I_file),@(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = imread(i_p.Results.cell_mask);
end 

%Pull out the parameters specfied on the command line
debug = i_p.Results.debug;
filt_size = i_p.Results.filt_size;
min_intensity = i_p.Results.min_intensity;
output_dir = i_p.Results.output_dir;

%read in and normalize the input focal adhesion image
focal_image  = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_filt = fspecial('average',filt_size);
high_passed_image = focal_image - imfilter(focal_image,I_filt,'same',mean(focal_image(:)));

adhesions = zeros(size(focal_image,1),size(focal_image,2));
adhesions(high_passed_image > min_intensity) = 1;
adhesions = imfill(adhesions,'holes');
if (exist('cell_mask','var'))
    adhesions = find_in_cell_ads(adhesions,cell_mask);
end

if (exist('cell_mask','var'))
    [B,F,T] = otsuThresholding(high_passed_image(cell_mask));
    adhesions_otsu = im2bw(high_passed_image,T);
    adhesions_otsu = find_in_cell_ads(adhesions_otsu,cell_mask);
    imwrite(adhesions_otsu,fullfile(output_dir, 'adhesions_otsu.png'));
end

if (exist('cell_mask','var'))
    adhesion_properties = collect_adhesion_properties(adhesions,focal_image,'cell_mask',cell_mask);
else
    adhesion_properties = collect_adhesion_properties(adhesions,focal_image);
end
    
%write the results to files
imwrite(adhesions,fullfile(output_dir, 'adhesions.png'));
write_adhesion_data(adhesion_properties,'out_dir',fullfile(output_dir,'raw_data'));

if (nargout > 0)
    varargout{1} = adhesions;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered_adhesions = find_in_cell_ads(ad,cm)
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
    
    labeled_ad = bwlabel(ad,4);
    
    filtered_adhesions = zeros(size(ad,1),size(ad,2));
    
    for i = 1:max(labeled_ad(:))
        this_adhesion = zeros(size(ad,1),size(ad,2));
        this_adhesion(labeled_ad == i) = 1;
        overlap = and(this_adhesion,cm);
        
        if (sum(overlap(:)))
            filtered_adhesions(overlap) = 1;
        end
    end
    
end

function adhesion_props = collect_adhesion_properties(ad_I,orig_I,varargin)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'COLLECT_ADHESION_PROPERTIES';

i_p.addRequired('ad_I',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('orig_I',@isnumeric);

i_p.addParamValue('cell_mask',0,@(x)isnumeric(x) || islogical(x));

i_p.parse(ad_I,orig_I,varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = i_p.Results.cell_mask;
end 

labeled_adhesions = bwlabel(ad_I,4);
adhesion_props = regionprops(labeled_adhesions,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Average_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Variance_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Max_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Min_adhesion_signal = zeros(1,max(labeled_adhesions(:)));

for i=1:max(labeled_adhesions(:))
    Average_adhesion_signal(i) = mean(orig_I(labeled_adhesions == i));
    Variance_adhesion_signal(i) = var(orig_I(labeled_adhesions == i));
    Max_adhesion_signal(i) = max(orig_I(labeled_adhesions == i));
    Min_adhesion_signal(i) = min(orig_I(labeled_adhesions == i));
end

adhesion_props(1).Average_adhesion_signal = Average_adhesion_signal;
adhesion_props(1).Average_pixel_adhesion_signal = Average_adhesion_signal*(2^12 - 1);
adhesion_props(1).Variance_adhesion_signal = Variance_adhesion_signal;
adhesion_props(1).Max_adhesion_signal = Max_adhesion_signal;
adhesion_props(1).Min_adhesion_signal = Min_adhesion_signal;

if (exist('cell_mask','var'))
    dists = bwdist(~cell_mask);
    cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
    cell_centroid = cell_centroid.Centroid;

    Centroid_dist_from_edge = zeros(1,max(labeled_adhesions(:)));
    Centroid_dist_from_center = zeros(1,max(labeled_adhesions(:)));
    Angle_to_center = zeros(1,max(labeled_adhesions(:)));

    for i=1:max(labeled_adhesions(:))
        Average_adhesion_signal(i) = mean(orig_I(labeled_adhesions == i));
        Variance_adhesion_signal(i) = var(orig_I(labeled_adhesions == i));

        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid;
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            Centroid_dist_from_edge(i) = dists(centroid_pos(2),centroid_pos(1));

            Centroid_dist_from_center(i) = sqrt((cell_centroid(1) - centroid_unrounded(1))^2 + (cell_centroid(2) - centroid_unrounded(2))^2);
            Angle_to_center(i) = acos((centroid_unrounded(2) - cell_centroid(2))/Centroid_dist_from_center(i));
            if (centroid_pos(2) - cell_centroid(2) < 0)
                Angle_to_center(i) = Angle_to_center(i) + pi;
            end
        end

    end

    adhesion_props(1).Centroid_dist_from_edge = Centroid_dist_from_edge;
    adhesion_props(1).Centroid_dist_from_center = Centroid_dist_from_center;
    adhesion_props(1).Angle_to_center = Angle_to_center;
    adhesion_props(1).Cell_size = sum(cell_mask(:));
end

end

end
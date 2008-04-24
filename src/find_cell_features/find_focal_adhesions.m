function [varargout] = find_focal_adhesions(I_file,cell_mask,varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,C,OPTIONS) locate the focal adhesions in image
%   file, 'I', and collects the properties of the focal adhesions using the
%   cell mask in file , 'C'
%
%   Options:
%
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

i_p.addRequired('cell_mask',@(x)exist(x,'file') == 2);

i_p.addParamValue('filt_size',23,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('min_intensity',0.1,@isnumeric);
i_p.addParamValue('output_dir',fileparts(I_file),@(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,cell_mask,varargin{:});

%Pull out the parameters specfied on the command line
debug = i_p.Results.debug;
filt_size = i_p.Results.filt_size;
min_intensity = i_p.Results.min_intensity;
output_dir = i_p.Results.output_dir;

%read in and normalize the input focal adhesion image
focal_image  = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%read in the cell mask image
cell_mask = imread(cell_mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_filt = fspecial('average',filt_size);
high_passed_image = focal_image - imfilter(focal_image,I_filt,'same',mean(focal_image(:)));

adhesions = zeros(size(focal_image,1),size(focal_image,2));
adhesions(find(high_passed_image > min_intensity)) = 1;
adhesions = imfill(adhesions,'holes');
adhesions = find_in_cell_ads(adhesions,cell_mask);

[B,F,T] = otsuThresholding(high_passed_image(find(cell_mask)));
adhesions_otsu = im2bw(high_passed_image,T);
adhesions_otsu = find_in_cell_ads(adhesions_otsu,cell_mask);

adhesion_properties = collect_adhesion_properties(adhesions,cell_mask,focal_image);

class_1_ind = find([adhesion_properties.Class] == 1);
class_2_ind = find([adhesion_properties.Class] == 2);
h = plot(adhesion_properties(1).Centroid_dist_from_center(class_1_ind),adhesion_properties(1).Centroid_dist_from_edge(class_1_ind),'ro');
hold on;
xlabel('Distance from Center (pixels)','FontSize',16);
ylabel('Distance from Edge (pixels)','FontSize',16);
plot(adhesion_properties(1).Centroid_dist_from_center(class_2_ind),adhesion_properties(1).Centroid_dist_from_edge(class_2_ind),'go');
plot(mean(adhesion_properties(1).Centroid_dist_from_center(class_2_ind)),mean(adhesion_properties(1).Centroid_dist_from_edge(class_2_ind)),'k+','MarkerSize',20,'LineWidth',4);
plot(mean(adhesion_properties(1).Centroid_dist_from_center(class_1_ind)),mean(adhesion_properties(1).Centroid_dist_from_edge(class_1_ind)),'k+','MarkerSize',20,'LineWidth',4);
hold off;
saveas(h,fullfile(output_dir, 'class_plot.pdf'));
saveas(h,fullfile(output_dir, 'class_plot.png'));
close all;

labeled_ad = bwlabel(adhesions,4);
class_highlight = focal_image;

for i = 1:max(labeled_ad(:))
    this_ad = zeros(size(labeled_ad,1),size(labeled_ad,2));
    this_ad(find(labeled_ad == i)) = 1;
    this_ad = bwperim(this_ad);
    class_highlight = create_highlighted_image(class_highlight,this_ad,'color',adhesion_properties(1).Class(i));
end
imwrite(class_highlight,fullfile(output_dir, 'ad_class.png'));

%write the results to files
imwrite(adhesions,fullfile(output_dir, 'adhesions.png'));
imwrite(adhesions_otsu,fullfile(output_dir, 'adhesions_otsu.png'));
write_adhesion_data(adhesion_properties,'out_dir',fullfile(output_dir,'raw_data'));

if (nargout > 0)
    varargout{1} = adhesions;
end

function filtered_adhesions = find_in_cell_ads(ad,cm)
    
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

end
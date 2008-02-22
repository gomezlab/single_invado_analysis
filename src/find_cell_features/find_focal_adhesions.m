function [varargout] = find_focal_adhesions(I_file,cell_mask)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,'out_dir',OUT_D,EXTRA_OPTIONS) finds the focal
%   adhesions in image 'I', which is a single image file, the resulting
%   binary segmented image is output in 'OUT_D' using the name
%   'focal_adhesions.png', the parameters in 'EXTRA_OPTIONS' will also be
%   used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_FOCAL_ADHESIONS';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('cell_mask',@(x)exist(x,'file') == 2);

i_p.parse(I_file,cell_mask);

%read in and normalize the input focal adhesion image
focal_image = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image   = double(focal_image)/scale_factor;

%check to see if a cell_mask parameter has been set, if it has, use that
%file, otherwise, search for the file in the output dir
cell_mask = imread(cell_mask);

output_dir = fileparts(I_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%focal_markers = find_focal_adhesion_markers(focal_image,cell_mask);
%watershed_edges = locate_watershed_edges(focal_image,focal_markers,cell_mask);

filt_size = 23;

I_filt = (1/filt_size^2)*ones(filt_size,filt_size);
high_passed_image = focal_image - conv2(focal_image,I_filt,'same');

min_filtered_pix_intensity = 0.1;

adhesions = zeros(size(focal_image));
adhesions(find(high_passed_image > min_filtered_pix_intensity)) = 1;

%adhesions = find_watershed_adhesions(focal_image,watershed_edges,cell_mask);
adhesion_properties = collect_adhesion_properties(adhesions,cell_mask,focal_image);

%write the results to files
%imwrite(focal_markers,fullfile(output_dir,'focal_markers.png'));
%imwrite(watershed_edges,fullfile(output_dir, 'watershed_edges.png'));
imwrite(adhesions,fullfile(output_dir, 'adhesions.png'));
write_adhesion_data(adhesion_properties,'out_dir',fullfile(output_dir,'raw_data'));

adhesion_props_filename = fullfile(output_dir, 'adhesion_props.mat');
adhesion_properties = adhesion_properties;
save(adhesion_props_filename, 'adhesion_properties');

if (nargout > 0)
    varargout{1} = adhesions;
end

end

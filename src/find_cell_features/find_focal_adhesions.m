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
%       -filt_size:
%       -min_intensity:
%       -output_dir:
%       -debug:

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

%Pull out the 
filt_size = i_p.Results.filt_size;
min_intensity = i_p.Results.min_intensity;
debug = i_p.Results.debug;
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

adhesion_properties = collect_adhesion_properties(adhesions,cell_mask,focal_image);

%write the results to files
imwrite(adhesions,fullfile(output_dir, 'adhesions.png'));
write_adhesion_data(adhesion_properties,'out_dir',fullfile(output_dir,'raw_data'));

if (nargout > 0)
    varargout{1} = adhesions;
end

end

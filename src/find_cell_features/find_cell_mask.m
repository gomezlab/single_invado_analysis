function [varargout] = find_cell_mask(I_file,out_file)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image 
%
%   create_cell_mask_image(I,OF) finds the cell mask using the image in
%   file 'I' and writes the binary cell mask to the output file 'OF'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('out_file',@(x)isempty(fileparts(x)) == 1 || exist(fileparts(x),'dir') == 7);

i_p.parse(I_file,out_file);

mask_image   = imread(I_file);
scale_factor = double(intmax(class(mask_image)));
mask_image   = double(mask_image)/scale_factor;

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sorted_mask_pixels = sort(mask_image(:));
% sorted_mask_pixels(1:0.05*round(length(sorted_mask_pixels))) = 0;

[heights, intensity] = hist(sorted_mask_pixels,1000);

smoothed_heights = smooth(heights,0.05,'loess');
[zmax,imax,zmin,imin]= extrema(smoothed_heights);

sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

min_indexes = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));

if (length(min_indexes) > 1)
    min_index = min_indexes(1)
    for i=2:length(min_indexes)
        if (imin(i) < min_index) 
            min_index = min_indexes(i);
        end
    end
else 
    min_index = min_indexes;
end

threshed_mask = im2bw(mask_image, intensity(imin(min_index)));

threshed_mask = clean_up_mask_image(threshed_mask);
imwrite(threshed_mask, out_file)

if (nargout >= 1)
    varargout{1} = threshed_mask;
end
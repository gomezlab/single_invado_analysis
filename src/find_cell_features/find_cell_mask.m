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
i_p.addRequired('out_file',@(x)exist(fileparts(x),'dir') == 7);

i_p.parse(I_file,out_file);

mask_image   = imread(I_file);
scale_factor = double(intmax(class(mask_image)));
mask_image   = double(mask_image)/scale_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_binary_image = bwperim(im2bw(mask_image,adaptive_thresh(mask_image,'upper_mean_weight',0.2)));
mask_binary_image = clean_up_mask_image(mask_binary_image);

imwrite(mask_binary_image,out_file);

if (nargout >= 1)
    varargout{1} = mask_binary_image;
end
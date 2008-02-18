function [varargout] = find_cell_mask(I_file,out_file)
%CREATE_CELL_MASK_IMAGE   Gather the cell mask from a fluorescence image
%
%   create_cell_mask_image(I) finds the cell mask using the image in file
%   'I', which contains only a single image
%
%   create_cell_mask_image(I,'out_dir',OD) finds the cell mask using the
%   image in file 'I', which contains only a single image and writes the
%   binary cell mask to the output directory 'OD'
%
%   create_cell_mask_image(I,'I_num',NUM) finds the cell mask using image
%   number 'NUM' in file 'I', which contains a multi-image stack
%
%   create_cell_mask_image(I,'I_num',NUM,'out_dir',OD) finds the cell mask
%   using image number 'NUM' in file 'I', which contains a multi-image
%   stack and writes the binary cell mask to the output directory 'OD'
%
%   Note: If an output parameter is requested the cell mask binary image is
%   returned.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('I_file',@(x)exist(x,'file')==2);

i_p.addRequired('out_file',pwd,@ischar);

i_p.parse(I_file,varargin{:});

%determine if the out_dir was set in the parameters, if set, fill in the
%output_dir field
if (not(any(strcmp('out_dir',i_p.UsingDefaults))))
    output_dir = i_p.Results.out_dir;
    if (not(exist(output_dir,'dir')))
        mkdir(output_dir);
    end
end



output_file = 'cell_mask.png';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_binary_image = bwperim(im2bw(mask_image,adaptive_thresh(mask_image,'upper_mean_weight',0.2)));
mask_binary_image = clean_up_mask_image(mask_binary_image);
mask_binary_image = imfill(mask_binary_image,'holes');
mask_binary_image = imdilate(mask_binary_image,strel('diamond',1));

if (exist('output_dir','var'))
    imwrite(mask_binary_image,fullfile(output_dir,output_file));
end

if (nargout >= 1)
    varargout{1} = mask_binary_image;
end

function [varargout] = create_cell_mask_image(varargin)
%CREATE_CELL_MASK_IMAGE   Gather the cell mask from a fluorescene image
%
%   create_cell_mask_image(I) extracts the cell mask from image 'I',
%   NOTE: the cell mask image will NOT be saved on the hard drive, only
%         returned to the calling function
%
%   create_cell_mask_image(I_FILE) extracts the cell mask from image
%   'I_FILE'
%   NOTE: that this file must only contain a single image
%   NOTE: the cell mask image will NOT be saved on the hard drive, only
%         returned to the calling function
%
%   create_cell_mask_image(I_FILE,I_NUM) extracts the cell mask from image
%   'I_FILE' with image number 'I_NUM'
%
%   create_cell_mask_image(I,OUT_DIR) extracts the cell mask from image 'I'
%   and writes the binary cell mask to the output directory 'OUT_DIR', with
%   the file name 'cell_mask.png'
%
%   create_cell_mask_image(I,OUT_DIR,OUT_FILE) extracts the cell mask from
%   image 'I' and writes the binary cell mask to the output directory
%   'OUT_DIR', with the file name 'OUT_FILE'
%
%   create_cell_mask_image(I_FILE,OUT_DIR) extracts the cell mask from file
%   'I_FILE' and writes the binary cell mask to the output directory
%   'OUT_DIR', with the file name 'cell_mask.png'
%
%   create_cell_mask_image(I_FILE,OUT_DIR,OUT_FILE) extracts the cell mask
%   from file 'I_FILE' and writes the binary cell mask to the output
%   directory 'OUT_DIR', with the file name 'OUT_FILE'
%
%   create_cell_mask_image(I_FILE,I_NUM,OUT_DIR) extracts the cell mask
%   from file 'I_FILE' with image number 'I_NUM' and writes the binary cell
%   mask to the output directory 'OUT_DIR', with the file name
%   'cell_mask.png'
%
%   create_cell_mask_image(I_FILE,I_NUM,OUT_DIR,OUT_FILE) extracts the cell mask
%   from file 'I_FILE' with image number 'I_NUM'
%   and writes the binary cell mask to the output directory 'OUT_DIR', with
%   the file name 'OUT_FILE'
%
%   If an output parameter is requested the cell mask binary image is
%   returned.

%%Setup Parameters and Variables
if (ischar(varargin{1}))
    first_input_type = exist(varargin{1});
else
    if (isnumeric(varargin{1}))
        first_input_type = 1;
    else
        error('ERROR: create_cell_mask_image - first parameter is not of the expected type');
    end
end

output_file = 'cell_mask.png';
write_output_files = true;

switch nargin
    case 1
        write_output_files = false;
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                mask_image = varargin{1};
            else
                error('ERROR: create_cell_mask_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            mask_image = imread(varargin{1});
        else
            error('ERROR: create_cell_mask_image - first parameter is not of the expected type');
        end
    case 2
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                mask_image = varargin{1};
                output_directory = varargin{2};
            else
                error('ERROR: create_cell_mask_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            mask_image = imread(varargin{1});
            if (isnumeric(varargin{2}))
                write_output_files = false;
                mask_image = imread(varargin{1},varargin{2});
            else
                output_directory = varargin{2};
            end
        else
            error('ERROR: create_cell_mask_image - first parameter is not of the expected type');
        end
    case 3
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                mask_image = varargin{1};
                output_directory = varargin{2};
                output_file = varargin{3};

            else
                error('ERROR: create_cell_mask_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            picture_info = imfinfo(varargin{1});
            if (size(picture_info,2) > 1)
                mask_image = imread(varargin{1},varargin{2});
                output_directory = varargin{3};
            else
                mask_image = imread(varargin{1});
                output_directory = varargin{2};
                output_file = varargin{3};
            end
        else
            error('ERROR: create_cell_mask_image - first parameter is not of the expected type');
        end
    case 4
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                mask_image = varargin{1};
                output_directory = varargin{2};
                output_file = varargin{3};
            else
                error('ERROR: create_cell_mask_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            picture_info = imfinfo(varargin{1});
            if (size(picture_info,2) > 1)
                mask_image = imread(varargin{1},varargin{2});
                output_directory = varargin{3};
                output_file = varargin{4};
            else
                mask_image = imread(varargin{1});
                output_directory = varargin{2};
                output_file = varargin{3};
            end
        else
            error('ERROR: create_cell_mask_image - first parameter is not of the expected type');
        end
    otherwise
        error('ERROR: create_cell_mask_image - unexpected number of arguments, expected 2-4');
end

if (min(mask_image(:)) ~= 0 && max(mask_image(:)) ~= 1)
    mask_image = normalize_grayscale_image(mask_image);
end

if (write_output_files)
    if (output_directory(end) ~= '/')
        output_directory = [output_directory, '/'];
    end
end

%%Main Program
mask_binary_image = bwperim(im2bw(mask_image,adaptive_thresh(mask_image,0.2)));
mask_binary_image = clean_up_mask_image(mask_binary_image);
mask_binary_image = imfill(mask_binary_image,'holes');
mask_binary_image = imdilate(mask_binary_image,strel('diamond',1));

if (write_output_files)
    if (not(exist(output_directory,'dir')))
        mkdir(output_directory);
    end

    imwrite(mask_binary_image,[output_directory,output_file]);
end

if (nargout >= 1)
    varargout{1} = mask_binary_image;
end
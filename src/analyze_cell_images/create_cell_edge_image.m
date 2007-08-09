function [varargout] = create_cell_edge_image(varargin)
%CREATE_CELL_EDGE_IMAGE   Gather the cell edge from a fluorescene image
%
%   create_cell_edge_image(I,OUT_DIR) extracts the cell edge from image 'I'
%   and writes the binary cell mask to the output directory 'OUT_DIR', with
%   the file name 'cell_mask.png'
%
%   create_cell_edge_image(I,OUT_DIR,OUT_FILE) extracts the cell edge from
%   image 'I' and writes the binary cell mask to the output directory
%   'OUT_DIR', with the file name 'OUT_FILE'
%
%   create_cell_edge_image(I_FILE,OUT_DIR) extracts the cell edge from file
%   'I_FILE' and writes the binary cell mask to the output directory
%   'OUT_DIR', with the file name 'cell_mask.png'
%
%   create_cell_edge_image(I_FILE,OUT_DIR,OUT_FILE) extracts the cell edge
%   from file 'I_FILE' and writes the binary cell mask to the output
%   directory 'OUT_DIR', with the file name 'OUT_FILE'
%
%   create_cell_edge_image(I_FILE,I_NUM,OUT_DIR) extracts the cell edge
%   from file 'I_FILE' with image number 'I_NUM' and writes the binary cell
%   mask to the output directory 'OUT_DIR', with the file name
%   'cell_mask.png'
%
%   create_cell_edge_image(I_FILE,I_NUM,OUT_DIR,OUT_FILE) extracts the cell edge
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
        error('ERROR: create_cell_edge_image - first parameter is not of the expected type');
    end
end

output_file = 'cell_mask.png';

switch nargin
    case 2
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                edge_image = varargin{1};
                output_directory = varargin{2};
            else
                error('ERROR: create_cell_edge_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            edge_image = imread(varargin{1});
            output_directory = varargin{2};
        else
            error('ERROR: create_cell_edge_image - first parameter is not of the expected type');
        end
    case 3
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                edge_image = varargin{1};
                output_directory = varargin{2};
                output_file = varargin{3};

            else
                error('ERROR: create_cell_edge_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            picture_info = imfinfo(varargin{1});
            if (size(picture_info,2) > 1)
                edge_image = imread(varargin{1},varargin{2});
                output_directory = varargin{3};
            else
                edge_image = imread(varargin{1});
                output_directory = varargin{2};
                output_file = varargin{3};
            end
        else
            error('ERROR: create_cell_edge_image - first parameter is not of the expected type');
        end
    case 4
        if (first_input_type == 1)
            if (isnumeric(varargin{1}))
                edge_image = varargin{1};
                output_directory = varargin{2};
                output_file = varargin{3};
            else
                error('ERROR: create_cell_edge_image - the first parameter is a matlab variable, but not a numeric matrix');
            end
        elseif (first_input_type == 2)
            picture_info = imfinfo(varargin{1});
            if (size(picture_info,2) > 1)
                edge_image = imread(varargin{1},varargin{2});
                output_directory = varargin{3};
                output_file = varargin{4};
            else
                edge_image = imread(varargin{1});
                output_directory = varargin{2};
                output_file = varargin{3};
            end
        else
            error('ERROR: create_cell_edge_image - first parameter is not of the expected type');
        end
    otherwise
        error('ERROR: create_cell_edge_image - unexpected number of arguments, expected 2-4');
end

if (min(edge_image(:)) ~= 0 && max(edge_image(:)) ~= 1)
    edge_image = normalize_grayscale_image(edge_image);
end

if (output_directory(end) ~= '/')
    output_directory = [output_directory, '/'];
end

%%Main Program
edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));
edge_binary_image = clean_up_edge_image(edge_binary_image);
edge_binary_image = imfill(edge_binary_image,'holes');
if (not(exist(output_directory,'dir')))
    mkdir(output_directory);
end

imwrite(edge_binary_image,[output_directory,output_file]);

if (nargout >= 1)
    varargout{1} = edge_binary_image;
end
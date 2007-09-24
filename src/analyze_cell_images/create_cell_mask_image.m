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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_file = 'cell_mask.png';
write_output_files = true;

if (isempty(varargin))
    error('ERROR: create_cell_mask_image - No parameters specified, see ''help create_cell_mask_image''');
elseif (isnumeric(varargin{1}))
    mask_image = varargin{1};
    if (length(varargin) == 1)
        write_output_files = false;
    elseif (length(varargin) == 2)
        if (not(exist(varargin{2},'dir')))
            error('ERROR: create_cell_mask_image - If the first parameter is a variable, expected second parameter to specify a directory');
        end
            output_directory = varargin{2};
    elseif (length(varargin) == 3)
        if (not(exist(varargin{2},'dir')))
            error('ERROR: create_cell_mask_image - If the first parameter is a variable, expected second parameter to specify a directory');            
        end
        output_directory = varargin{2};
        output_file = varargin{3};
    else
        error('ERROR: create_cell_mask_image - If the first parameter is a numeric, expected 1, 2 or 3 parameters');
    end
elseif (exist(varargin{1},'file'))
    if (length(varargin) == 1)
        write_output_files = false;
        mask_image = imread(varargin{1});
    else
        if (isnumeric(varargin{2}))
            if (length(varargin) == 2)
                write_output_files = false;
                mask_image = imread(varargin{1},varargin{2});
            elseif (length(varargin) == 3)
                mask_image = imread(varargin{1},varargin{2});
                if (exist(varargin{3},'dir'))
                    output_directory = varargin{3};
                else
                    error('ERROR: create_cell_mask_image - First parameter was a file, second a number, third parameter must be a directory');
                end
            elseif (length(varargin) == 4)
                mask_image = imread(varargin{1},varargin{2});

                if (not(exist(varargin{3},'dir')))
                    error('ERROR: create_cell_mask_image - First parameter was a file, second a number, third parameter must be a directory');
                end
                if (not(ischar(varargin{4})))
                    error('ERROR: create_cell_mask_image - First parameter was a file, second a number, third parameter a directory, fourth should be a the desired filename');
                end

                output_directory = varargin{3};
                output_file = varargin{4};
            else
                error('ERROR: create_cell_mask_image - expected parameter counts of 2, 3 or 4 when first input is a file and second a number');
            end
        elseif (exist(varargin{2},'dir'))
            if (length(varargin) == 2)
                mask_image = imread(varargin{1});
                output_directory = varargin{2};
            elseif (length(varargin) == 3)
                mask_image = imread(varargin{1});
                output_directory = varargin{2};

                if (not(ischar(varargin{3})))
                    error('ERROR: create_cell_mask_image - First parameter was a file, second a directory, thrird parameter should be the desired filename');
                end

                output_file = varargin{3};

            end
        else
            error('ERROR: create_cell_mask_image - First parameter was a file, second should be either a directory or a the image number');
        end
    end
else
    error('ERROR: create_cell_mask_image - Problem with the parameters, see ''help create_cell_mask_image''');
end

if (min(mask_image(:)) ~= 0 && max(mask_image(:)) ~= 1)
    mask_image = normalize_grayscale_image(mask_image);
end

if (write_output_files)
    if (output_directory(end) ~= '/')
        output_directory = [output_directory, '/'];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_binary_image = bwperim(im2bw(mask_image,adaptive_thresh(mask_image,0.2)));
mask_binary_image = clean_up_mask_image(mask_binary_image);
mask_binary_image = imfill(mask_binary_image,'holes');
mask_binary_image = imdilate(mask_binary_image,strel('diamond',1));

if (write_output_files)
    if (not(exist(output_directory,'dir')))
        mkdir(output_directory);
    end

    imwrite(mask_binary_image,fullfile(output_directory,output_file));
end

if (nargout >= 1)
    varargout{1} = mask_binary_image;
end

function [varargout]=find_focal_adhesions(varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,OUT_DIR,EXTRA_OPTIONS) finds the focal adhesions
%   in image 'I', which can be either a single image file or a matlab
%   variable, the resulting binary segmented image is output in 'OUT_DIR'
%   using the name 'focal_adhesions.png', the parameters in 'EXTRA_OPTIONS'
%   will also be used
%
%   find_focal_adhesions(I,EXTRA_OPTIONS) finds the focal adhesions
%   in image 'I', which can be either a single image file or a matlab
%   variable, the location of the cell mask must be specified in
%   'EXTRA_OPTIONS'
%
%   find_focal_adhesions(I,I_NUM,OUT_DIR,EXTRA_OPTIONS) finds the focal
%   adhesions in file 'I' image number 'I_NUM', the resulting binary
%   segmented image is output in 'OUT_DIR' using the name
%   'focal_adhesions.png', the parameters in 'EXTRA_OPTIONS' will also be
%   used
%
%   find_focal_adhesions(I,I_NUM,EXTRA_OPTIONS) finds the focal
%   adhesions in file 'I' image number 'I_NUM', the parameters in
%   'EXTRA_OPTIONS' will also be used, the location of the cell mask must
%   be specified in 'EXTRA_OPTIONS'
%
%   NOTE: The cell mask, which is the binary image indicating the inside
%   versus the outside of the cell must be present in the output directory,
%   with the file name 'cell_mask.png', otherwise, it must be specified in
%   the extra options section.
%
%   NOTE: If a return variable is requested, the binary segmented focal
%   adhesions will be returned
%
%   Extra Options:
%       -'debug' - if present, debug mode on, yielding dianostic messages,
%        short form 'd'
%       -'cell_mask' - the next command line parameter must either be a
%        matlab variable with the cell mask or a file where the mask can be
%        found, short form 'c'
%
%
%   Example Commands:
%
%       find_focal_adhesions('/path/to/image.png','/output/dir/')
%           -assumes that the cell mask is available in the output dir
%           '/output/dir'
%
%       find_focal_adhesions(I,'cell_mask','/path/to/mask.png')
%           -finds the focal adhesions in the matlab variable 'I', using
%           the cell_mask found in '/path/to/mask.png'
%
%       find_focal_adhesions(I,'c','/path/to/mask.png')
%            -equivalent to the last example, only using the short form of
%            the 'cell_mask' extra option

%%Setup variables and parse command line

image_data = struct;
image_data.debug = 0;

%check to see if the command has enough to parameters to be valid
if (length(varargin) < 3)
    error('ERROR: find_focal_adhesions - not enough command line options, see ''help find_focal_adhesions''');
end

%now parse the non-'extra options' part of the command line
if (exist(varargin{1},'file'))
    if (isnumeric(varargin{2}))
        image_data.original_image = imread(varargin{1},varargin{2});
        if (exist(varargin{3},'dir'))
            image_data.output_dir = varargin{3};
        end
    elseif (exist(varargin{2},'dir'))
        image_data.original_image = imread(varargin{1});
        image_data.output_dir = varargin{2};
    end
elseif (exist(varargin{1},'var'))
    image_data.original_image = varargin{1};
    if (exist(varargin{2},'dir'))
        image_data.output_dir = varargin{2};
    end
end

%now parse the 'extra options'
extra_options_parse = struct;
extra_options_parse.cell_mask_position = 0;

for i = 1:length(varargin)
    if (ischar(varargin{i}))
        if (strcmpi(varargin{i},'cell_mask') || strcmpi(varargin{i},'c'))
            if (length(varargin) <= i)
                error('ERROR: find_focal_adhesions - must specify the cell mask location directly after the command line option')
            end
            extra_options_parse.cell_mask = varargin{i + 1};
            extra_options_parse.cell_mask_position = i + 1;
        end
        if (strcmpi(varargin{i},'debug') || strcmpi(varargin{i},'d'))
            image_data.debug = 1;
        end
    end
end

if (extra_options_parse.cell_mask_position)
    if (not(islogical(extra_options_parse.cell_mask) || isnumeric(extra_options_parse.cell_mask) || exist(extra_options_parse.cell_mask,'file')))
        error('ERROR: find_focal_adhesions - the cell mask must be either a matlab variable, of logical or numeric type, or a file')
    else
        if (islogical(extra_options_parse.cell_mask) || isnumeric(extra_options_parse.cell_mask))
            image_data.cell_mask = extra_options_parse.cell_mask;
        elseif (exist(extra_options_parse.cell_mask,'file'))
            image_data.cell_mask = imread(extra_options_parse.cell_mask);
        end
    end
end

%now do some final logical checks to make sure all the variables needed are
%available

if (not(isfield(image_data,'output_dir')) && not(isfield(image_data,'cell_mask')))
    error('ERROR: find_focal_adhesions - without an output dir the cell mask must be specified in the extra options');
end

if (not(isfield(image_data,'cell_mask')))
    poss_cell_mask_file = fullfile(image_data.output_dir, 'cell_mask.png');
    if (not(exist(poss_cell_mask_file,'file')))
        error('ERROR: find_focal_adhesions - cannot find the cell mask in either the extra options or in the output directory under the file ''cell_mask.png''')
    else
        image_data.cell_mask = imread(poss_cell_mask_file);
    end
end

%now normalize the input focal adhesion image
image_data.original_image = normalize_grayscale_image(image_data.original_image);

%%Main Program

image_data.focal_markers = find_focal_adhesion_markers(image_data);


if (nargout > 0)
    varargout{1} = image_data;
end

end
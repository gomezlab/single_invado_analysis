function [varargout] = find_focal_adhesions_val(I_file,varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,OUT_DIR,EXTRA_OPTIONS) finds the focal adhesions
%   in image 'I', which is a single image file, the resulting binary
%   segmented image is output in 'OUT_DIR' using the name
%   'focal_adhesions.png', the parameters in 'EXTRA_OPTIONS' will also be
%   used
%
%   find_focal_adhesions(I,EXTRA_OPTIONS) finds the focal adhesions
%   in image 'I', which is a single image file, the location of the cell
%   mask must be specified in 'EXTRA_OPTIONS' 
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
%   NOTE: If a return variable is requested, the struct holding the results
%   of many of the steps in the algorithm is returned
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_data = struct;
image_data.debug = 0;

i_p = inputParser;
i_p.FunctionName = 'FIND_FOCAL_ADHESIONS';

i_p.addRequired('I_file',@(x)exist(x,'file')==2);

i_p.addParamValue('out_dir',pwd,@ischar);
i_p.addParamValue('I_num',1,@(x)isnumeric(x) && x>0);
i_p.addParamValue('debug',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));
i_p.addParamValue('d',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));
i_p.addParamValue('cell_mask',0,@(x)exist(x,'file')==2);
i_p.addParamValue('c',0,@(x)exist(x,'file')==2);

i_p.parse(I_file,varargin{:});

%Start filling out the image_data struct with the parameters extracted
%above

%Pull out the original image file and the data directory
image_data.original_image_file = i_p.Results.I_file;
[pathstr, name, ext, versn] = fileparts(image_data.original_image_file); 
image_data.data_folder = pathstr;

%Determine if the image file specified has more than one image embedded, if
%so, make sure there is a 'I_num' parameter
if (size(imfinfo(image_data.original_image_file),2) > 1)
    if (any(strcmp('I_num',i_p.UsingDefaults)))
        error('ERROR: FIND_FOCAL_ADHESIONS - Image file specified has more than one image embedded, must specify an ''I_num'' parameter');
    end
    image_data.I_num = i_p.Results.I_num;
end

%Set the debug flag on if either 'd' or 'debug' flag has been set to true
image_data.debug = i_p.Results.d || i_p.Results.debug;

%determine if the out_dir was set in the parameters, if set, fill in the
%output_dir field
if (not(any(strcmp('out_dir',i_p.UsingDefaults))))
    image_data.output_dir = i_p.Results.out_dir;
    if (not(exist(image_data.output_dir,'dir')))
        mkdir(image_data.output_dir);
    end
end

%check to see if a min/max value file has been written, if not collect the
%min/max values and write them to file
if (not(exist(fullfile(image_data.data_folder,'min_max_vals.txt'),'file')))
    [min,max] = find_extr_values(image_data.original_file);
    extr_vals = [min,max];
    save(fullfile(image_data.data_folder,'min_max_vals.txt'),'extr_vals','-ASCII');
else
    extr_vals = load(fullfile(image_data.data_folder,'min_max_vals.txt'));
end
image_data.extr_vals = extr_vals;

%normalize the input focal adhesion image
if (isfield(image_data,'I_num'))
    image_data.original_image = imread(image_data.original_image_file,image_data.I_num);
else
    image_data.original_image = imread(image_data.original_image_file);
end
image_data.original_image = normalize_grayscale_image(image_data.original_image,image_data.extr_vals(1),image_data.extr_vals(2));

%check to see if a cell_mask parameter has been set, if it has, use that
%file, otherwise, search for the file in the output dir
if (not(any(strcmp('c',i_p.UsingDefaults))) || not(any(strcmp('cell_mask',i_p.UsingDefaults))))
    if (not(strcmp('cell_mask',i_p.UsingDefaults)))
        image_data.cell_mask_file = i_p.Results.cell_mask;
    else
        image_data.cell_mask_file = i_p.Results.c;
    end
elseif (isfield(image_data,'output_dir'))
    if (exist(fullfile(image_data.output_dir,'cell_mask.png'),'file'))
        image_data.cell_mask_file = fullfile(image_data.output_dir,'cell_mask.png');
    else
        error('ERROR: find_focal_adhesions - could not locate cell mask file, specify location with ''cell_mask'' flag');
    end
else
    error('ERROR: find_focal_adhesions - could not locate cell mask file, specify location with ''cell_mask'' flag');
end
image_data.cell_mask = imread(image_data.cell_mask_file);

if (image_data.debug)
    disp(' ')
    disp 'List of all arguments:'
    disp(i_p.Results)

    disp 'Parameters Using Defaults:'
    disp(i_p.UsingDefaults)

    disp 'Contents of image_data:'
    disp(image_data)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_data.focal_markers = find_focal_adhesion_markers(image_data.original_image,image_data.cell_mask);
image_data.watershed_edges = locate_watershed_edges(image_data);
image_data.adhesions = find_watershed_adhesions(image_data);
image_data.adhesion_properties = collect_adhesion_properties(image_data);

if (isfield(image_data,'output_dir'))
    imwrite(image_data.original_image,fullfile(image_data.output_dir, 'focal_image.png'));
    imwrite(image_data.focal_markers,fullfile(image_data.output_dir,'focal_markers.png'));
    imwrite(image_data.watershed_edges,fullfile(image_data.output_dir, 'watershed_edges.png'));
    imwrite(image_data.adhesions,fullfile(image_data.output_dir, 'adhesions.png'));
    write_adhesion_data(image_data.adhesion_properties,'dir',image_data.output_dir);
    
    adhesion_props_filename = fullfile(image_data.output_dir, 'adhesion_props.mat');
    adhesion_properties = image_data.adhesion_properties;
    save(adhesion_props_filename, 'adhesion_properties');
end

if (nargout > 0)
    varargout{1} = image_data;
end

end
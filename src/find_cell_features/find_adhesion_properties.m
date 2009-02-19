function find_adhesion_properties(focal_file,adhesions_file,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps
%
%   find_adhesion_properties(ff,af,OPTIONS) determines the quantitative
%   properites of the adhesions identified in the file 'af', using
%   information in the focal image file 'ff', the properties are written
%   to a set of csv files located a folder named 'raw_data' in the same
%   directory as the focal image file 'ff'
%
%   Options:
%
%       -cell_mask: file which contains the cell mask, defaults to not
%        present
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_ADHESION_PROPERTIES';

i_p.addRequired('focal_file',@(x)exist(x,'file') == 2);
i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);

i_p.parse(focal_file, adhesions_file);

i_p.addParamValue('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addOptional('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addOptional('protrusion_file',0,@(x)exist(x,'file') == 2);
i_p.addOptional('i_num',0,@(x)isnumeric(x) && x >= 1);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image  = imread(focal_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%read in the labeled adhesions
adhesions = imread(adhesions_file);

%check if protrusion_file is specified, read it in if i_num is also
%specified
if (isempty(strmatch('protrusion_file',i_p.UsingDefaults)))
    if (not(strmatch('i_num',i_p.UsingDefaults)))
        warning('AdhesionProps:dataloading','Protrusion file specified, but don''t know which image number to use for this image, not including protrusion data.');
    else
        load(i_p.Results.protrusion_file);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    if (exist('protrusion','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'protrusion_data',protrusion,'i_num',i_p.Results.i_num,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'debug',i_p.Results.debug);
    end
else
    if (exist('protrusion_matrix','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'protrusion_data',protrusion,'i_num',i_p.Results.i_num,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'debug',i_p.Results.debug);
    end
end
if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));
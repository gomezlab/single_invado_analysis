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

i_p.addOptional('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addOptional('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image  = imread(focal_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

adhesions = imread(adhesions_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    adhesion_properties = collect_adhesion_properties(adhesions,focal_image,'cell_mask',cell_mask);
else
    adhesion_properties = collect_adhesion_properties(adhesions,focal_image);
end
if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));
function find_cell_mask_properties(current_dir,first_dir,final_dir,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;

i_p = inputParser;

i_p.addRequired('current_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('first_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('final_dir',@(x)exist(x,'dir') == 7);

i_p.parse(current_dir,first_dir, final_dir);

i_p.addParamValue('output_dir',current_dir,@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(current_dir, first_dir, final_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_data = struct;

%read in and normalize the input focal adhesion image
current_data.gel_image  = imread(fullfile(current_dir, filenames.gel_filename));
scale_factor = double(intmax(class(current_data.gel_image)));
current_data.gel_image  = double(current_data.gel_image)/scale_factor;

%read in the intensity correction coefficient
current_data.intensity_correction = csvread(fullfile(current_dir, filenames.intensity_correction_filename));

%read in the cell mask file
current_data.cell_mask = logical(imread(fullfile(current_dir, filenames.cell_mask_filename)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the final directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_data = struct;

%read in and normalize the input focal adhesion image
final_data.gel_image  = imread(fullfile(final_dir, filenames.gel_filename));
scale_factor = double(intmax(class(final_data.gel_image)));
final_data.gel_image  = double(final_data.gel_image)/scale_factor;

%read in the labeled adhesions
final_data.binary_shift = logical(imread(fullfile(final_dir, filenames.binary_shift_filename)));

%read in the intensity correction coefficient
final_data.intensity_correction = csvread(fullfile(final_dir, filenames.intensity_correction_filename));

%read in the cell mask file
final_data.cell_mask = logical(imread(fullfile(final_dir, filenames.cell_mask_filename)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the first directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_data = struct;

%read in and normalize the input focal adhesion image
first_data.gel_image  = imread(fullfile(first_dir, filenames.gel_filename));
scale_factor = double(intmax(class(first_data.gel_image)));
first_data.gel_image  = double(first_data.gel_image)/scale_factor;

first_data.binary_shift = logical(imread(fullfile(first_dir, filenames.binary_shift_filename)));
first_data.intensity_correction = csvread(fullfile(first_dir, filenames.intensity_correction_filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_data.adhesions = bwlabel(current_data.cell_mask,8);
final_data.adhesions = bwlabel(final_data.cell_mask,8);

adhesion_properties = collect_adhesion_properties(current_data, first_data, final_data,'debug',i_p.Results.debug);

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));

profile off;
if (i_p.Results.debug), profile viewer; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = collect_adhesion_properties(c_d,first_d,f_d,varargin)
% COLLECT_ADHESION_PROPERTIES    using the identified adhesions, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                adhesions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_ADHESION_PROPERTIES';

i_p.addRequired('c_d',@isstruct);
i_p.addRequired('first_d',@isstruct);
i_p.addRequired('f_d',@isstruct);

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(c_d,first_d,f_d,varargin{:});

adhesion_props = regionprops(c_d.adhesions,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(c_d.adhesions(:))    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = c_d.adhesions;
    this_ad(c_d.adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    
%     differences = c_d.gel_image(this_ad)*c_d.intensity_correction - first_d.gel_image(this_ad)*first_d.intensity_correction;
    differences = c_d.gel_image(this_ad) - first_d.gel_image(this_ad);
    [h,p] = ttest(differences);
    
    adhesion_props(i).Cell_gel_diff_p_val = p;
    adhesion_props(i).Cell_gel_diff = mean(differences);
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(c_d.adhesions(:)))]); end
end
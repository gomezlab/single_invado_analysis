function find_adhesion_properties(current_dir,first_dir,final_dir,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_PUNCTA_PROPERTIES';

i_p.addRequired('current_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('first_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('final_dir',@(x)exist(x,'dir') == 7);

i_p.parse(current_dir,first_dir, final_dir);

i_p.addParamValue('output_dir',current_dir,@ischar);
i_p.addParamValue('adhesions_filename','puncta_labeled.png',@ischar);
i_p.addParamValue('puncta_filename','registered_focal_image.png',@ischar);
i_p.addParamValue('gel_filename','registered_gel.png',@ischar);
i_p.addParamValue('binary_shift_filename','binary_shift.png',@ischar);
i_p.addParamValue('cell_mask_filename','cell_mask.png',@ischar);
i_p.addParamValue('intensity_correction_file','intensity_correction.csv',@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(current_dir, first_dir, final_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_data = struct;

%read in and normalize the input focal adhesion image
current_data.puncta_image  = imread(fullfile(current_dir, i_p.Results.puncta_filename));
scale_factor = double(intmax(class(current_data.puncta_image)));
current_data.puncta_image  = double(current_data.puncta_image)/scale_factor;

%read in and normalize the input focal adhesion image
current_data.gel_image  = imread(fullfile(current_dir, i_p.Results.gel_filename));
scale_factor = double(intmax(class(current_data.gel_image)));
current_data.gel_image  = double(current_data.gel_image)/scale_factor;

%read in the labeled adhesions
current_data.adhesions = imread(fullfile(current_dir, i_p.Results.adhesions_filename));

%read in the labeled adhesions
current_data.binary_shift = logical(imread(fullfile(current_dir, i_p.Results.binary_shift_filename)));

%read in the intensity correction coefficient
current_data.intensity_correction = csvread(fullfile(current_dir, i_p.Results.intensity_correction_file));

%read in the cell mask file if defined
if(exist(fullfile(current_dir, i_p.Results.cell_mask_filename), 'file'))
    current_data.cell_mask = logical(imread(fullfile(current_dir, i_p.Results.cell_mask_filename)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the final directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_data = struct;

%read in and normalize the input focal adhesion image
final_data.puncta_image  = imread(fullfile(final_dir, i_p.Results.puncta_filename));
scale_factor = double(intmax(class(final_data.puncta_image)));
final_data.puncta_image  = double(final_data.puncta_image)/scale_factor;

%read in and normalize the input focal adhesion image
final_data.gel_image  = imread(fullfile(final_dir, i_p.Results.gel_filename));
scale_factor = double(intmax(class(final_data.gel_image)));
final_data.gel_image  = double(final_data.gel_image)/scale_factor;

%read in the labeled adhesions
final_data.adhesions = imread(fullfile(final_dir, i_p.Results.adhesions_filename));

%read in the labeled adhesions
final_data.binary_shift = logical(imread(fullfile(final_dir, i_p.Results.binary_shift_filename)));

%read in the intensity correction coefficient
final_data.intensity_correction = csvread(fullfile(final_dir, i_p.Results.intensity_correction_file));

%read in the cell mask file if defined
if(exist(fullfile(final_dir, i_p.Results.cell_mask_filename), 'file'))
    final_data.cell_mask = logical(imread(fullfile(final_dir, i_p.Results.cell_mask_filename)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the final directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_data = struct;

%read in and normalize the input focal adhesion image
first_data.gel_image  = imread(fullfile(first_dir, i_p.Results.gel_filename));
scale_factor = double(intmax(class(first_data.gel_image)));
first_data.gel_image  = double(first_data.gel_image)/scale_factor;

%read in the labeled adhesions
first_data.adhesions = imread(fullfile(first_dir, i_p.Results.adhesions_filename));

%read in the labeled adhesions
first_data.binary_shift = logical(imread(fullfile(first_dir, i_p.Results.binary_shift_filename)));

%read in the intensity correction coefficient
first_data.intensity_correction = csvread(fullfile(first_dir, i_p.Results.intensity_correction_file));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adhesion_properties = collect_adhesion_properties(current_data, first_data, final_data,'debug',i_p.Results.debug);

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));

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
    adhesion_props(i).Average_adhesion_signal = mean(c_d.puncta_image(c_d.adhesions == i));
    adhesion_props(i).Variance_adhesion_signal = var(c_d.puncta_image(c_d.adhesions == i));
    adhesion_props(i).Max_adhesion_signal = max(c_d.puncta_image(c_d.adhesions == i));
    adhesion_props(i).Min_adhesion_signal = min(c_d.puncta_image(c_d.adhesions == i));
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = c_d.adhesions;
    this_ad(c_d.adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    
    current_diffs = collect_local_diff_properties(c_d,this_ad);
    
    adhesion_props(i).Local_gel_diff = current_diffs.Local_gel_diff;
    adhesion_props(i).Local_gel_diff_corr = current_diffs.Local_gel_diff_corr;
    adhesion_props(i).Global_gel_diff = current_diffs.Global_gel_diff;
    adhesion_props(i).Large_local_gel_diff = current_diffs.Large_local_gel_diff;
    adhesion_props(i).Large_local_gel_diff_corr = current_diffs.Large_local_gel_diff_corr;
    
    first_diffs = collect_local_diff_properties(first_d,this_ad);
    adhesion_props(i).First_local_gel_diff = first_diffs.Local_gel_diff;
    
    final_diffs = collect_local_diff_properties(f_d,this_ad);
    adhesion_props(i).End_local_gel_diff = final_diffs.Local_gel_diff;
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(c_d.adhesions(:)))]); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(c_d, 'cell_mask'))
    [dists, indicies] = bwdist(~c_d.cell_mask);
    dists(not(c_d.binary_shift)) = NaN;
    
    min_row = find(sum(c_d.binary_shift,2),1,'first');
    max_row = find(sum(c_d.binary_shift,2),1,'last');
    min_col = find(sum(c_d.binary_shift),1,'first');
    max_col = find(sum(c_d.binary_shift),1,'last');
    
    only_reg_mask = c_d.cell_mask(min_row:max_row, min_col:max_col);
    assert(size(only_reg_mask,1)*size(only_reg_mask,2) == sum(sum(c_d.binary_shift)));
    
    %Now we search for the pixels which are closest to an edge of the cell
    %mask that is also touching the edge of image. We want to find these
    %pixels because the true closest cell edge may be off the microscope
    %field of view. To be safe, we will set those
    %distance-to-nearest-cell-edge values to NaN.
    black_border_mask = only_reg_mask;
    black_border_mask(1,:) = 0; black_border_mask(end,:) = 0;
    black_border_mask(:,1) = 0; black_border_mask(:,end) = 0;
    
    bb_dists_temp = bwdist(~black_border_mask);
    
    bb_dists = zeros(size(c_d.cell_mask));
    bb_dists(min_row:max_row, min_col:max_col) = bb_dists_temp;
    dists(not(c_d.binary_shift)) = NaN;
    
    dists(bb_dists < dists) = NaN;
    for i=1:max(c_d.adhesions(:))
        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid;
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
        end
    end
end
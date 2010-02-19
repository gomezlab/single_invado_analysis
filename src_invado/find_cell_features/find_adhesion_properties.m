function find_adhesion_properties(current_dir,final_dir,varargin)
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
i_p.addRequired('final_dir',@(x)exist(x,'dir') == 7);

i_p.parse(current_dir, final_dir);

i_p.addParamValue('output_dir',current_dir,@ischar);
i_p.addParamValue('adhesions_filename','puncta_labeled.png',@ischar);
i_p.addParamValue('puncta_filename','registered_focal_image.png',@ischar);
i_p.addParamValue('gel_filename','registered_gel.png',@ischar);
i_p.addParamValue('binary_shift_filename','binary_shift.png',@ischar);
i_p.addParamValue('cell_mask_filename','cell_mask.png',@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(current_dir, final_dir,varargin{:});

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

%read in the cell mask file if defined
if(exist(fullfile(final_dir, i_p.Results.cell_mask_filename), 'file'))
    final_data.cell_mask = logical(imread(fullfile(final_dir, i_p.Results.cell_mask_filename)));
end

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adhesion_properties = collect_adhesion_properties(current_data, final_data,'debug',i_p.Results.debug);

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = collect_adhesion_properties(c_d,f_d,varargin)
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
i_p.addRequired('f_d',@isstruct);

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(c_d,f_d,varargin{:});

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
    background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
    %we don't want to include any areas that have been identified as other
    %objects
    background_region = and(background_region,not(c_d.adhesions));
    %also exclude areas outside the registered image
    background_region = logical(background_region .* c_d.binary_shift);
    assert(sum(sum(background_region)) > 0)
    
    adhesion_props(i).Local_gel_diff = mean(c_d.gel_image(this_ad)) - mean(c_d.gel_image(background_region));
    adhesion_props(i).Global_gel_diff = mean(c_d.gel_image(this_ad)) - mean(c_d.gel_image(c_d.binary_shift));
    
    adhesion_props(i).Background_adhesion_signal = mean(c_d.puncta_image(background_region));
    adhesion_props(i).Background_area = sum(background_region(:));
    adhesion_props(i).Background_corrected_signal = adhesion_props(i).Average_adhesion_signal - adhesion_props(i).Background_adhesion_signal;
    
    large_background_region = logical(imdilate(this_ad,strel('disk',10,0)));
    %we don't want to include any areas that have been identified as other
    %objects
    large_background_region = and(large_background_region,not(c_d.adhesions));
    %also exclude areas outside the registered image
    large_background_region = logical(large_background_region .* c_d.binary_shift);
    assert(sum(sum(large_background_region)) > 0)
    
    adhesion_props(i).Large_local_gel_diff = mean(c_d.gel_image(this_ad)) - mean(c_d.gel_image(large_background_region));
    
    final_background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
    final_background_region = and(final_background_region,not(f_d.adhesions));
    final_background_region = logical(final_background_region .* f_d.binary_shift);
    if (sum(sum(final_background_region)) > 0) 
        adhesion_props(i).End_local_gel_diff = mean(f_d.gel_image(this_ad)) - mean(f_d.gel_image(final_background_region));
    else 
        adhesion_props(i).End_local_gel_diff = NaN;
    end
    
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

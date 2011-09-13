function find_puncta_properties(current_dir,first_dir,final_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

i_p = inputParser;
i_p.FunctionName = 'FIND_PUNCTA_PROPERTIES';

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
% Pull in data from directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_data = read_in_file_set(current_dir,filenames);
final_data = read_in_file_set(final_dir,filenames);
first_data = read_in_file_set(first_dir,filenames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adhesion_properties = collect_adhesion_properties(current_data, first_data, final_data,'debug',i_p.Results.debug);

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));

toc;

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

adhesion_props = regionprops(c_d.objects,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(c_d.objects(:))
    adhesion_props(i).Average_puncta_signal = mean(c_d.puncta_image(c_d.objects == i));
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = c_d.objects;
    this_ad(c_d.objects ~= i) = 0;
    this_ad = logical(this_ad);
    
    current_diffs = collect_local_diff_properties(c_d,this_ad);
    
    adhesion_props(i).Local_gel_diff = current_diffs.Local_gel_diff;
    adhesion_props(i).Global_gel_diff = current_diffs.Global_gel_diff;
    adhesion_props(i).Large_local_gel_diff = current_diffs.Large_local_gel_diff;
    
    first_diffs = collect_local_diff_properties(first_d,this_ad);
    adhesion_props(i).First_local_gel_diff = first_diffs.Local_gel_diff;
    
    final_diffs = collect_local_diff_properties(f_d,this_ad);
    adhesion_props(i).End_local_gel_diff = final_diffs.Local_gel_diff;
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(c_d.objects(:)))]); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(c_d, 'cell_mask'))
    [dists, indicies] = bwdist(~c_d.cell_mask); %#ok<NASGU>
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
    for i=1:max(c_d.objects(:))
        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid; %#ok<NASGU>
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
        end
    end
end
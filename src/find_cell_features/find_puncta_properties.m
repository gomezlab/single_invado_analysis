function find_puncta_properties(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

first_data = read_in_file_set(fullfile(base_dir,image_dirs(1).name),filenames);
final_data = read_in_file_set(fullfile(base_dir,image_dirs(end).name),filenames);

for i_num = 1:size(image_dirs,1)
    current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    object_properties = collect_object_properties(current_data, first_data, final_data,'debug',i_p.Results.debug);
    
    %write the results to files
    write_object_data(object_properties,'out_dir',fullfile(base_dir,image_dirs(i_num).name,'raw_data'));
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function object_props = collect_object_properties(c_d,first_d,f_d,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('c_d',@isstruct);
i_p.addRequired('first_d',@isstruct);
i_p.addRequired('f_d',@isstruct);

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(c_d,first_d,f_d,varargin{:});

object_props = regionprops(c_d.objects,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(c_d.objects(:))
    object_props(i).Average_puncta_signal = mean(c_d.puncta_image(c_d.objects == i));
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = c_d.objects;
    this_ad(c_d.objects ~= i) = 0;
    this_ad = logical(this_ad);
    
    current_diffs = collect_local_diff_properties(c_d,this_ad);
    
    object_props(i).Local_gel_diff = current_diffs.Local_gel_diff;
    object_props(i).Global_gel_diff = current_diffs.Global_gel_diff;
    object_props(i).Large_local_gel_diff = current_diffs.Large_local_gel_diff;
    
    first_diffs = collect_local_diff_properties(first_d,this_ad);
    object_props(i).First_local_gel_diff = first_diffs.Local_gel_diff;
    
    final_diffs = collect_local_diff_properties(f_d,this_ad);
    object_props(i).End_local_gel_diff = final_diffs.Local_gel_diff;
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(c_d.objects(:)))]); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isfield(c_d, 'cell_mask'))
    [dists, ~] = bwdist(~c_d.cell_mask);
    
    %Now we search for the pixels which are closest to an edge of the cell
    %mask that is also touching the edge of image. We want to find these
    %pixels because the true closest cell edge may be off the microscope
    %field of view. To be safe, we will set those
    %distance-to-nearest-cell-edge values to NaN.
    black_border_mask = ones(size(c_d.cell_mask,1),size(c_d.cell_mask,2));
    black_border_mask(1,:) = 0; black_border_mask(end,:) = 0;
    black_border_mask(:,1) = 0; black_border_mask(:,end) = 0;
    
    bb_dists_temp = bwdist(~black_border_mask);
        
    dists(bb_dists_temp < dists) = NaN;
    for i=1:max(c_d.objects(:))
        centroid_pos = round(object_props(i).Centroid);
        centroid_unrounded = object_props(i).Centroid; %#ok<NASGU>
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_object_properties - centroid not found');
            object_props(i).Centroid_dist_from_edge = NaN;
        else
            object_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
        end
    end
end

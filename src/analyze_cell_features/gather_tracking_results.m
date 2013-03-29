function gather_tracking_results(field_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('field_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('pixel_size',NaN,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(field_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));
% addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(field_dir,'individual_pictures');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

tracking_file = fullfile(base_dir, image_dirs(1).name,filenames.tracking);
%check for the existance of a tracking file, if absent, there weren't any
%cells in this field, return from the function
if (not(exist(tracking_file,'file')))
    disp('No tracking matrix found, assuming no cells in field');
    return;
else
    tracking_mat = csvread(tracking_file) + 1;
end

load(fullfile(base_dir,image_dirs(1).name,filenames.cell_props));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the types of data in set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_types = fieldnames(all_cell_props{1}); %#ok<USENS>
for i_num = 2:length(all_cell_props)
	%need this loop in case the first image set was empty
	if (isempty(data_types)) 
		data_types = fieldnames(all_cell_props{1});
	else
		continue;
	end
end

data_to_exclude = {'MeanIntensity','StdIntensity','Centroid','Overlap_area',...
    'Overlap_percent','Eccentricity'};

output_dir = fullfile(base_dir,image_dirs(1).name,filenames.lineage_dir);
if (not(exist(output_dir,'dir'))), mkdir(output_dir); end
    
for i = 1:length(data_types)
    if (any(strcmp(data_types{i},data_to_exclude)))
        continue;
    end
    
    this_data = scan_through_property(tracking_mat,all_cell_props,data_types{i});
    
    if (strcmp(data_types{i},'Area') && not(isnan(i_p.Results.pixel_size)))
        this_data = this_data * i_p.Results.pixel_size^2;
    end
    
    dlmwrite(fullfile(output_dir,[data_types{i}, '.csv']), this_data);
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results_mat = scan_through_property(tracking_mat,all_cell_props,prop_name)

results_mat = NaN(size(tracking_mat));

for tracking_num = 1:size(tracking_mat,1)
    for i_num = 1:size(tracking_mat,2)
        cell_num = tracking_mat(tracking_num,i_num);
        if (cell_num < 1)
            continue;
        end
        
        results_mat(tracking_num,i_num) = all_cell_props{i_num}(cell_num).(prop_name);
    end
end

function [results_mat_1, results_mat_2] = scan_through_bi_property(tracking_mat,all_cell_props,prop_name)

results_mat_1 = NaN(size(tracking_mat));
results_mat_2 = NaN(size(tracking_mat));

for tracking_num = 1:size(tracking_mat,1)
    for i_num = 1:size(tracking_mat,2)
        cell_num = tracking_mat(tracking_num,i_num);
        if (cell_num < 1)
            continue;
        end
        
        results_mat_1(tracking_num,i_num) = all_cell_props{i_num}(cell_num).(prop_name)(1);
        results_mat_2(tracking_num,i_num) = all_cell_props{i_num}(cell_num).(prop_name)(2);
        
    end
end

function cell_speed = determine_cell_speed(centroid_x,centroid_y)
    
cell_speed = NaN(size(centroid_x));

for cell_num = 1:size(centroid_x,1)
    for i_num = 2:size(centroid_x,2)
        if (not(isnan(centroid_x(cell_num,i_num - 1))) && ...
            not(isnan(centroid_x(cell_num,i_num))))
            
            speed = sqrt((centroid_x(cell_num,i_num) - centroid_x(cell_num, i_num - 1))^2 + ...
                (centroid_y(cell_num,i_num) - centroid_y(cell_num, i_num - 1))^2);
            
            cell_speed(cell_num,i_num) = speed;
        end
    end
end

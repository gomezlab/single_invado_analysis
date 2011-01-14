function create_invader_visualization(field_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('field_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(field_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

active_degrade = csvread(fullfile(field_dir,'adhesion_props','active_degrade.csv'));

%add one because perl indexes from 0
tracking_matrix = csvread(fullfile(field_dir,'tracking_matrices','tracking_seq.csv')) + 1;

base_dir = fullfile(field_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

for i_num = 1:size(image_dirs)
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    labeled_cells_perim = zeros(size(current_data.labeled_cells));
    for cell_num=1:max(current_data.labeled_cells(:))
        temp = zeros(size(current_data.labeled_cells));
        temp(current_data.labeled_cells == cell_num) = 1;
        temp = bwperim(temp);
        labeled_cells_perim(temp) = cell_num;
    end
    
    degrade_marked = zeros(size(current_data.labeled_cells));
    for cell_num=1:max(labeled_cells_perim(:))
        this_tracking_row = tracking_matrix(:,i_num) == cell_num;
        assert(sum(this_tracking_row) == 1);
        
        degrade_status = active_degrade(this_tracking_row,i_num);
        
        %use 1 for non-degraders, 2 for degraders
        if (degrade_status == 0)
            degrade_marked(labeled_cells_perim == cell_num) = 1;
        elseif (degrade_status == 1)
            degrade_marked(labeled_cells_perim == cell_num) = 2;
        else
            disp('Found unexpected degrade status');
        end
    end
    
    centroid = csvread(fullfile(base_dir,image_dirs(i_num).name,'raw_data','Centroid.csv'));
    gel_diff = csvread(fullfile(base_dir,image_dirs(i_num).name,'raw_data','Cell_gel_diff.csv'));
    gel_diff_percent = csvread(fullfile(base_dir,image_dirs(i_num).name,'raw_data','Cell_gel_diff_percent.csv'));
    
    norm_gel = (current_data.gel_no_norm - current_data.gel_range(1))/range(current_data.gel_range);
    
    c_map = [[1,0,0];[0,1,0]];
    
    degrade_highlights = create_highlighted_image(norm_gel,degrade_marked,'color_map',c_map);
    
    imshow(degrade_highlights);
    for cell_num = 1:length(gel_diff)
        text(centroid(cell_num,1), centroid(cell_num,2), ... 
            fprintf('%d\n%d',gel_diff(cell_num),gel_diff_percent(cell_num)));
    end
    
    1;
end
toc;

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

data_sets_to_read = {'Cell_gel_diff','Cell_gel_diff_median','Cell_gel_diff_p_val','Area' ...
    'Cell_gel_diff_total','Centroid_x','Centroid_y'};
raw_data = struct();

for i = 1:length(data_sets_to_read)
    data_dir = fullfile(base_dir, image_dirs(1).name,filenames.lineage_dir);
    raw_data.(data_sets_to_read{i}) = csvread(fullfile(data_dir,[data_sets_to_read{i}, '.csv']));
end

tracking_mat = csvread(fullfile(base_dir, image_dirs(1).name,filenames.tracking))+1;
%occastionally, there will be fields were no cells were detected, in that
%case, the tracking matrix is an empty file, but we still need a matrix to
%index into, so make an empty tracking matrix with a line of zeros
if (isempty(tracking_mat))
    tracking_mat = zeros(1,size(image_dirs,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building Visualization Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 1:size(image_dirs,1)
    output_file = fullfile(base_dir,image_dirs(i_num).name, filenames.invader_vis);
    if (not(i_p.Results.debug) && exist(output_file,'file'))
        continue;
    end
    
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    degrade_marked = zeros(size(current_data.labeled_cells));
    for cell_num=1:max(current_data.labeled_cells_perim(:))
        this_tracking_row = tracking_matrix(:,i_num) == cell_num;
        assert(sum(this_tracking_row) == 1);
        
        degrade_status = active_degrade(this_tracking_row,i_num);
        
        %use 1 for non-degraders, 2 for degraders
        if (degrade_status == 0)
            degrade_marked(current_data.labeled_cells_perim == cell_num) = 1;
        elseif (degrade_status == 1)
            degrade_marked(current_data.labeled_cells_perim == cell_num) = 2;
        else
            disp('Found unexpected degrade status');
        end
    end

    c_map = [[1,0,0];[0,1,0]];
    
    degrade_highlights = create_highlighted_image(current_data.gel_image_norm,degrade_marked,'color_map',c_map);
    imwrite(degrade_highlights,output_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding Degradation Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convert_cmd = 'convert ';
convert_avail = not(system('which convert'));
if (not(convert_avail))
    if (exist('convert','file'))
        convert_cmd = './convert ';
    else
        exit();
    end
end
for i_num = 1:size(image_dirs,1)
    output_file = fullfile(base_dir,image_dirs(i_num).name, filenames.invader_vis);
    
    if (exist(output_file,'file'))
        img_size = size(imread(output_file));
        
        tracking_col = tracking_mat(:,i_num);
        
        filtered_data = struct();
        for i = 1:length(data_sets_to_read)
            temp = NaN(sum(tracking_col > 0),1);
            for j=1:length(tracking_col)
                if (tracking_col(j) < 1)
                    continue;
                end
                temp(tracking_col(j)) = raw_data.(data_sets_to_read{i})(j, i_num);
            end
            filtered_data.(data_sets_to_read{i}) = temp;
            
        end
        
        filtered_data.Centroid_x(filtered_data.Centroid_x > 0.9*img_size(2)) = 0.9*img_size(2);
        filtered_data.Centroid_y(filtered_data.Centroid_y > 0.9*img_size(1)) = 0.9*img_size(1);
        filtered_data.Centroid_y(filtered_data.Centroid_y < 0.1*img_size(1)) = 0.1*img_size(1);
        centroid = [filtered_data.Centroid_x,filtered_data.Centroid_y];
        
        gel_diff = filtered_data.Cell_gel_diff;
        area = filtered_data.Area;
        gel_diff_total = filtered_data.Cell_gel_diff_total;
        gel_diff_median = filtered_data.Cell_gel_diff_median;
        
        all_annotate = '';
        for cell_num = 1:length(gel_diff)
            pos_str = [' +',num2str(centroid(cell_num,1)),'+',num2str(centroid(cell_num,2))];
            label_str = [' "',sprintf('%.0f',area(cell_num)), '-test\n', ...
                sprintf('%.2f',gel_diff(cell_num)), '\n', ...
                sprintf('%.0f',gel_diff_total(cell_num)), '\n', ...
                sprintf('%.2f',gel_diff_median(cell_num)),'" '];
            all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
        end
        command_str = [convert_cmd, output_file, ' -font VeraBd.ttf -pointsize 16 -fill ''rgba(0,0,0,1)''', ...
            all_annotate, ' ', output_file, ';'];
        
        system(command_str);
    end
end

toc;

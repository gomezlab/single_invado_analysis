function build_single_cell_montage(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputparser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'error: expected the third string to be image set one') %#ok<st2nm>

image_dirs = image_dirs(3:end);

filenames = add_filenames_to_struct(struct());

tracking_seq = csvread(fullfile(base_dir,image_dirs(1).name,filenames.tracking)) + 1;
if (isempty(tracking_seq))
    tracking_seq = zeros(1,max_image_num);
end

for track_id = 1:size(tracking_seq,1)
    track_row = tracking_seq(track_id,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine extent of cell coverage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_coverage = [];
    for i_num = 1:length(track_row)
        if (track_row(i_num) <= 0)
            continue;
        end
        current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
        this_cell = current_data.labeled_cells == track_row(i_num);
        if (size(cell_coverage,1) == 0)
            cell_coverage = this_cell;
        else
            cell_coverage = cell_coverage | this_cell;
        end
    end
    
    row_min_max = [find(sum(cell_coverage,2),1,'first') - 10, ...
        find(sum(cell_coverage,2),1,'last') + 10];
    col_min_max = [find(sum(cell_coverage),1,'first') - 10,...
        find(sum(cell_coverage),1,'last') + 10];
    
    if (row_min_max(1) < 0), row_min_max(1) = 1; end
    if (col_min_max(1) < 0), col_min_max(1) = 1; end
    
    if (row_min_max(2) > size(cell_coverage,1)), row_min_max(2) = size(cell_coverage,1); end
    if (col_min_max(2) > size(cell_coverage,2)), col_min_max(2) = size(cell_coverage,2); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build visualization frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_frames = cell(0);
    puncta_frames = cell(0);
    for i_num = 1:length(track_row)
        if (track_row(i_num) <= 0)
            continue;
        end
        current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
        current_data = trim_all_images(current_data,row_min_max,col_min_max);        
        
        this_cell = current_data.labeled_cells == track_row(i_num);
        this_cell_perim = current_data.labeled_cells_perim == track_row(i_num);
        not_this_cell_perim = current_data.labeled_cells_perim ~= track_row(i_num) & ...
            current_data.labeled_cells_perim > 0;
        
        thick_perim = thicken_perimeter(this_cell_perim,this_cell);
        
        gel_frame = create_highlighted_image(current_data.gel_image_norm,thick_perim,'color_map',[1,0,0]);
        gel_frame = create_highlighted_image(gel_frame,not_this_cell_perim,'color_map',[0,0,0.5]);
        gel_frames{length(gel_frames) + 1} = gel_frame;
        
        puncta_frame = create_highlighted_image(current_data.puncta_image_norm,thick_perim,'color_map',[1,0,0]);
        puncta_frame = create_highlighted_image(puncta_frame,not_this_cell_perim,'color_map',[0,0,0.5]);
        puncta_frames{length(puncta_frames) + 1} = puncta_frame;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create final montage image and output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_montage = create_montage_image_set(gel_frames,'num_cols',length(gel_frames));
    puncta_montage = create_montage_image_set(puncta_frames,'num_cols',length(puncta_frames));
    spacer = 0.5*ones(1,size(gel_montage,2),3);
    
    full_montage = cat(1,puncta_montage,spacer,gel_montage);
    
    image_num = sprintf('%03d',track_id);
    output_file = fullfile(base_dir,image_dirs(1).name,filenames.single_cell_dir,[image_num,'.png']);
    
    if (not(exist(fileparts(output_file),'dir')))
        mkdir(fileparts(output_file))
    end
    
    imwrite(full_montage,output_file);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function images = trim_all_images(images,row_min_max,col_min_max)
image_size = size(images.gel_image);

names = fieldnames(images);
for i = 1:size(names)
    this_field = names{i};
    this_size = size(images.(this_field));
    if (length(this_size) == length(image_size) && ...
        all(this_size == image_size))
        
        1;
        temp = images.(this_field);
        temp = temp(row_min_max(1):row_min_max(2),...
            col_min_max(1):col_min_max(2));
        images.(this_field) = temp;
        1;
    end
end



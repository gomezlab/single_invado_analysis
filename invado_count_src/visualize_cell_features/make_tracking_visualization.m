function make_tracking_visualization(exp_dir,varargin)
%MAKE_MOVIE_FRAMES    Builds movie frames with the adhesions highlighted in
%                     using various conventions

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

filenames = add_filenames_to_struct(struct());

tracking_seq = csvread(fullfile(base_dir,image_dirs(1).name,filenames.tracking)) + 1;
if (isempty(tracking_seq))
    tracking_seq = zeros(1,max_image_num);
end

lineage_cmap = jet(size(tracking_seq,1));

for i_num = 1:length(image_dirs)
    current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);

    %thicken the cell perimeter, easier for visualization
    labeled_cells_perim_thick = current_data.labeled_cells_perim;
    three_disk = strel('disk',3,0);
    for cell_num = 1:max(labeled_cells_perim_thick(:))
        this_cell_perim = current_data.labeled_cells_perim == cell_num;
        this_cell = current_data.labeled_cells == cell_num;
        
        this_cell_perim = imdilate(this_cell_perim,three_disk);
        this_cell_perim = this_cell_perim & this_cell;
        
        labeled_cells_perim_thick(this_cell_perim) = cell_num;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    live_rows = tracking_seq(:,i_num) > 0;
    ad_nums = tracking_seq(live_rows,i_num);
    
    %Build the unique lineage highlighted image
    this_cmap(ad_nums,:) = lineage_cmap(live_rows,:); %#ok<AGROW>    
    highlighted_puncta = create_highlighted_image(current_data.puncta_image_norm, ...
        labeled_cells_perim_thick,'color_map',this_cmap);
    
    imwrite(highlighted_puncta,fullfile(base_dir,image_dirs(i_num).name,filenames.tracking_vis))
end
toc;
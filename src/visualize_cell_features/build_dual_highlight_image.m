function build_dual_highlight_image(exp_dir,varargin)
%BUILD_DUAL_HIGHLIGHT_IMAGES    Builds highlighted images showing the
%                               invadopodia (green), non-invadopodia (red),
%                               other puncta (blue) and the cell edge
%                               (purple), placing the image for each frame
%                               of the movie in the individual image
%                               folders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
i_p = inputParser;
i_p.FunctionName = 'BUILD_DUAL_HIGHLIGHT_IMAGES';

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('image_sets',NaN,@iscell);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'))
filenames = add_filenames_to_struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all the image directories
individual_images_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(individual_images_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invado_data_file = fullfile(individual_images_dir,image_dirs(1).name,filenames.invado_data);
try
    invado_data = csvread(invado_data_file,1,0);
catch err
    fprintf('Problem with invado data file: %s\n',err.message);
    invado_data = zeros(0,3);
end

not_invado_data_file = fullfile(individual_images_dir,image_dirs(1).name,filenames.not_invado_data);
try
    not_invado_data = csvread(not_invado_data_file,1,0);
catch err
    fprintf('Problem with not invado data file: %s\n',err.message);
    not_invado_data = zeros(0,3);
end

%after loading the tracking sequence filter to only include those puncta
%included in the invadopodia list, remember the list is formated so the first
%column contains the lineage number, with the first lineage as 1, so no need to
%translate
tracking_seq = load(fullfile(individual_images_dir,image_dirs(1).name,filenames.tracking_matrix)) + 1;
invado_tracking_seq = tracking_seq(invado_data(:,1),:);
not_invado_tracking_seq = tracking_seq(not_invado_data(:,1),:);

only_birth_point = zeros(size(tracking_seq,1),size(tracking_seq,2));
for puncta_num = 1:size(only_birth_point,1)
    birth_time = find(tracking_seq(puncta_num,:) > 0,1,'first');
    only_birth_point(puncta_num,birth_time) = tracking_seq(puncta_num,birth_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Reading, If Not All Ready Defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (not(any(strcmp(i_p.UsingDefaults,'image_sets'))))
    image_sets = i_p.Results.image_sets;
else
    image_sets = cell(size(image_dirs,1),1);
    for i = 1:size(image_dirs,1)
        image_sets{i} = read_in_file_set(fullfile(individual_images_dir,image_dirs(i).name),filenames);
    end
    toc(start_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_folder = fullfile(exp_dir,'visualizations','invado_and_not');
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

output_folder_side = fullfile(exp_dir,'visualizations','side_by_side');
if (not(exist(output_folder_side,'dir')))
    mkdir(output_folder_side);
end

output_folder_overlap = fullfile(exp_dir,'visualizations','overlap');
if (not(exist(output_folder_overlap,'dir')))
    mkdir(output_folder_overlap);
end

for i_num = 1:length(image_dirs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the object label perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    invado_nums = invado_tracking_seq(:,i_num);
    invado_nums = invado_nums(invado_nums > 0);
    
    not_invado_nums = not_invado_tracking_seq(:,i_num);
    not_invado_nums = not_invado_nums(not_invado_nums > 0);
    
    neither_nums = 1:max(image_sets{i_num}.objects_perim(:));
    neither_nums = setdiff(neither_nums,invado_nums);
    neither_nums = setdiff(neither_nums,not_invado_nums);
    
    assert(isempty(intersect(invado_nums,neither_nums)))
    assert(isempty(intersect(not_invado_nums,neither_nums)))
    assert(isempty(intersect(invado_nums,not_invado_nums)))
    
    puncta_label_perim_invado = ismember(image_sets{i_num}.objects_perim,invado_nums);
    puncta_label_perim_not_invado = ismember(image_sets{i_num}.objects_perim,not_invado_nums);
    puncta_label_perim_neither = ismember(image_sets{i_num}.objects_perim,neither_nums);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    images_to_highlight = {image_sets{i_num}.gel_image_norm, image_sets{i_num}.puncta_image_norm};
    original_images = images_to_highlight;
    spacer = ones(size(images_to_highlight{1},1),1,3);
    spacer_gray = ones(size(images_to_highlight{1},1),1);
    cell_perim = bwperim(image_sets{i_num}.cell_mask);
    
    for j=1:length(images_to_highlight)
        %cell edge highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},cell_perim,'color_map',[150/255,46/255,166/255]);
        
        %puncta that aren't classified
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_neither,'color_map',[0,0,1]);
        
        %not-invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_not_invado,'color_map',[1,0,0]);
        
        %invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_invado,'color_map',[0,1,0]);
    end
    
    output_image = [images_to_highlight{1}, spacer, images_to_highlight{2}];
    invado_and_not_output_image = fullfile(output_folder,[sprintf('%04d',i_num),'.png']);
    imwrite(output_image, invado_and_not_output_image);
    
    output_image = [original_images{1}, spacer_gray, original_images{2}];
    imwrite(output_image, fullfile(output_folder_side,[sprintf('%04d',i_num),'.png']));
    
    %Overlayed dual color image of puncta and ECM
    output_image = cat(3,original_images{1},...
        original_images{2},...
        zeros(size(original_images{1},1),size(original_images{1},2)));
    imwrite(output_image, fullfile(output_folder_overlap,[sprintf('%04d',i_num),'.png']));
    
    all_annotate = '';
    props = regionprops(image_sets{i_num}.objects);
    for puncta_num = [invado_nums;not_invado_nums]'
        lineage_num = find(only_birth_point(:,i_num) == puncta_num);
        if (isempty(lineage_num)), continue; end
        
        pos_str = [' +',num2str(props(puncta_num).Centroid(1)),'+',num2str(props(puncta_num).Centroid(2))];
        all_annotate = [all_annotate, ' -annotate', pos_str, ' ', num2str(lineage_num)]; %#ok<AGROW>
    end
    command_str = ['convert ', invado_and_not_output_image, ' -undercolor ''rgba(1,1,1,0.25)'' -font VeraBd.ttf -pointsize 16 -fill ''rgba(1,1,1,0.75)''', ...
        all_annotate, ' ', invado_and_not_output_image, ';'];
    system(command_str);
end

toc(start_time);
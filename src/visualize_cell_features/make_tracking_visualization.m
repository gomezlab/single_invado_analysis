function make_tracking_visualization(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('image_file','puncta',@ischar);
i_p.addParamValue('image_min_max_file','puncta_range',@ischar);
i_p.addParamValue('out_folder','tracking',@ischar);

i_p.addParamValue('pixel_size',0,@(x)x == 1 || x == 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('image_sets',NaN,@iscell);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));

filenames = add_filenames_to_struct(struct());

image_padding_min = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
individual_images_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(individual_images_dir);
image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Reading, If Not All Ready Defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (not(any(strcmp(i_p.UsingDefaults,'image_sets'))))
    image_sets = i_p.Results.image_sets;
else
    image_sets = cell(size(image_dirs,1),1);
    for i = 1:size(image_dirs,1)
        image_sets{i} = read_in_file_set(fullfile(individual_images_dir,image_dirs(i).name),filenames);
        if (mod(i,10) == 0)
            disp(['Finished Reading ', num2str(i), '/',num2str(size(image_dirs,1))]);
        end
    end
    toc(start_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in binary images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_binary = [];
for i_num = 1:length(image_dirs)
    this_binary = imread(fullfile(individual_images_dir,image_dirs(i_num).name,filenames.objects_binary));
    if (any(size(all_binary) == 0))
        all_binary = zeros(size(this_binary));
    end
    all_binary = all_binary | this_binary;
    all_binary = all_binary | image_sets{i_num}.cell_mask;
end
col_bounds = find(sum(all_binary));
col_bounds = [col_bounds(1) - image_padding_min,col_bounds(end) + image_padding_min];
row_bounds = find(sum(all_binary,2));
row_bounds = [row_bounds(1) - image_padding_min,row_bounds(end) + image_padding_min];

if (col_bounds(1) < 1), col_bounds(1) = 1; end
if (row_bounds(1) < 1), row_bounds(1) = 1; end
if (col_bounds(2) > size(all_binary,2)), col_bounds(2) = size(all_binary,2); end
if (row_bounds(2) > size(all_binary,1)), row_bounds(2) = size(all_binary,1); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign Each Object a Unique Color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracking_seq = csvread(fullfile(individual_images_dir,image_dirs(1).name,filenames.tracking)) + 1;

max_live_objects = max(sum(tracking_seq > 0));

lineage_cmap = jet(max_live_objects);
lineage_to_cmap = zeros(size(tracking_seq,1),1);
edge_cmap = jet(size(tracking_seq,2));


for i_num = 1:length(image_dirs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the matrices translating number to colormap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(tracking_seq,1)
        %if the object idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_num) <= 0), continue; end
        
        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_seq(:,i_num) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);
            
            taken_nums = zeros(1,max_live_objects);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);
            
            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
            catch map_error %#ok<NASGU>
                assert(~any(taken_dists == max(taken_dists)), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end
    
    %Make sure all the live objects have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_num) > 0) > 0), 'Error in assigning unique color codes');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ad_nums_lineage_order = tracking_seq(tracking_seq(:,i_num) > 0,i_num);
    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i_num) > 0);
    %     this_cmap = zeros(max(image_sets{i_num}.objects_perim(:)),3);
    this_cmap(ad_nums_lineage_order,:) = lineage_cmap(cmap_nums,:); %#ok<AGROW>
    highlighted_all = create_highlighted_image(image_sets{i_num}.puncta_image_norm,image_sets{i_num}.objects_perim,'color_map',this_cmap);
    
    cell_edge = bwperim(image_sets{i_num}.cell_mask);
    highlighted_all = create_highlighted_image(highlighted_all,cell_edge,'color_map',edge_cmap(i_num,:));
    
    %Bound the images according to the bounding box found towards the top
    %of the program
    image_sets{i_num}.puncta_image_norm = image_sets{i_num}.puncta_image_norm(row_bounds(1):row_bounds(2), col_bounds(1):col_bounds(2));
    highlighted_all = highlighted_all(row_bounds(1):row_bounds(2), col_bounds(1):col_bounds(2), 1:3);
    
    spacer = 0.5*ones(size(image_sets{i_num}.puncta_image_norm,1),1,3);
    
    composite_image = [cat(3,image_sets{i_num}.puncta_image_norm,...
                             image_sets{i_num}.puncta_image_norm,...
                             image_sets{i_num}.puncta_image_norm),...
                       spacer,...
                       highlighted_all];
    
    %Add scale bars if the pixel size is available
    if (exist('pixel_size','var'))
        composite_image = draw_scale_bar(composite_image,pixel_size);
    end
    
    out_folder = fullfile(exp_dir,'visualizations',i_p.Results.out_folder);
    if (not(exist(out_folder,'dir')))
        mkdir(out_folder);
    end
    
    out_file = fullfile(out_folder,sprintf('%05d.png',i_num));
    imwrite(composite_image,out_file);
    
    if (i_p.Results.debug)
        disp(['Done with ', num2str(i_num), '/',num2str(length(image_dirs))]);
    end
end

toc(start_time);
function make_tracking_visualization(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('image_file','puncta',@ischar);
i_p.addParamValue('image_min_max_file','puncta_range',@ischar);
i_p.addParamValue('out_folder','tracking',@ischar);

i_p.addParamValue('pixel_size',0,@(x)x == 1 || x == 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));

filenames = add_filenames_to_struct(struct());

image_padding_min = 10;

image_file = i_p.Results.image_file;
image_min_max_file = i_p.Results.image_min_max_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
individual_images_dir = fullfile(exp_dir,'individual_pictures');
image_folders = dir(individual_images_dir);
image_folders = image_folders(3:end);

%checking for presence of first file for visualization, if it isn't
%present, quit program
test_file = fullfile(individual_images_dir,image_folders(1).name,filenames.(image_file));
if (not(exist(test_file,'file')))
    exit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_binary = [];
for i_num = 1:length(image_folders)
    this_binary = imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.objects_binary));
    if (any(size(all_binary) == 0))
        all_binary = zeros(size(this_binary));
    end
    all_binary = all_binary | this_binary;
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
% Assign Each Adhesion a Unique Color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The tracking file won't exist in cases where a single image has been
%sent into the system. We still want to see the multi-color tracking-like
%visualization though, so build a dummy tracking_seq variable
if (exist(fullfile(exp_dir,filenames.tracking),'file'))
    tracking_seq = csvread(fullfile(exp_dir,filenames.tracking)) + 1;
else
    the_adhesions = imread(fullfile(individual_images_dir,image_folders(1).name,filenames.adhesions));
    tracking_seq = (1:max(the_adhesions(:)))';
end

max_live_adhesions = max(sum(tracking_seq > 0));

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);
edge_cmap = jet(size(tracking_seq,2));

image_set_range = csvread(fullfile(individual_images_dir,image_folders(1).name,filenames.(image_min_max_file)));
for i_num = 1:length(image_folders)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather and scale the input adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_i = double(imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.(image_file))));
    orig_i = (orig_i - image_set_range(1))/(image_set_range(2) - image_set_range(1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label_perim = imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.objects_perim));
    
    cell_mask_file = fullfile(individual_images_dir,image_folders(i_num).name,filenames.cell_mask);
    if (exist(cell_mask_file,'file'))
        cell_edge = bwperim(imread(cell_mask_file));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the matrices translating number to colormap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(tracking_seq,1)
        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_num) <= 0), continue; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_seq(:,i_num) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);

            taken_nums = zeros(1,max_live_adhesions);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);

            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
            catch map_error %#ok<NASGU>
                assert(~any(taken_dists == max(taken_dists)), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end

    %Make sure all the live adhesions have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_num) > 0) > 0), 'Error in assigning unique color codes');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ad_nums_lineage_order = tracking_seq(tracking_seq(:,i_num) > 0,i_num);
    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i_num) > 0);
%     this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums_lineage_order,:) = lineage_cmap(cmap_nums,:); %#ok<AGROW>
    highlighted_all = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);

    if (exist('cell_mask','var'))
        highlighted_all = create_highlighted_image(highlighted_all,cell_edge,'color_map',edge_cmap(i_num,:));
    end
    
    %Bound the images according to the bounding box found towards the top
    %of the program
    orig_i = orig_i(row_bounds(1):row_bounds(2), col_bounds(1):col_bounds(2));
    highlighted_all = highlighted_all(row_bounds(1):row_bounds(2), col_bounds(1):col_bounds(2), 1:3);

    spacer = 0.5*ones(size(orig_i,1),1,3);
    
    composite_image = [cat(3,orig_i,orig_i,orig_i),spacer,highlighted_all];

    %Add scale bars if the pixel size is available
    if (exist('pixel_size','var') && not(i_p.Results.no_scale_bar))
        composite_image = draw_scale_bar(composite_image,pixel_size);
    end
    
    out_folder = fullfile(exp_dir,'visualizations',i_p.Results.out_folder);
    if (not(exist(out_folder,'dir')))
        mkdir(out_folder);
    end
    
    out_file = fullfile(out_folder,sprintf('%05d.png',i_num));
    imwrite(composite_image,out_file);
    
    if (i_p.Results.debug)
        disp(['Done with ', num2str(i_num), '/',num2str(length(image_folders))]);
    end
end

toc;
profile off;
if (i_p.Results.debug), profile viewer; end

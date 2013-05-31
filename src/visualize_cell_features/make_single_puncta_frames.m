function make_single_puncta_frames(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('puncta_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('image_padding',15,@isnumeric);
i_p.addParamValue('image_sets',NaN,@iscell);
i_p.addParamValue('pixel_size',0,@(x)isnumeric(x));
i_p.addParamValue('sample_number',0,@(x)isnumeric(x));

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end
    toc(start_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking Sequence Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracking_seq = load(fullfile(individual_images_dir,image_dirs(1).name,filenames.tracking_matrix)) + 1;

invado_data_file = fullfile(individual_images_dir,image_dirs(1).name,filenames.invado_data);
not_invado_data_file = fullfile(individual_images_dir,image_dirs(1).name,filenames.not_invado_data);

tracking_temp = zeros(size(tracking_seq,1),size(tracking_seq,2));
try
    invado_data = csvread(invado_data_file,1,0);
    invado_nums = invado_data(:,1);
    %column one contains all the puncta lineage nums
    tracking_temp(invado_nums,:) = tracking_seq(invado_nums,:);
catch err
    fprintf('Problem with invado data file: %s\n',err.message);
end

try
    not_invado_data = csvread(not_invado_data_file,1,0);
    not_invado_nums = not_invado_data(:,1);
    %column one contains all the puncta lineage nums
    tracking_temp(not_invado_nums,:) = tracking_seq(not_invado_nums,:);
catch err
    fprintf('Problem with not invado data file: %s\n',err.message);
end
tracking_seq = tracking_temp;

if (exist('invado_nums','var') && exist('not_invado_nums','var'))
    assert(isempty(intersect(invado_nums,not_invado_nums)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Bounding Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The bounding matrix will hold the proper ranges to cut out only the
%relavent area for each object's life cycle. The matrix is initialized
%with values that will always be replaced and the data is structured like
%this:
%   Column 1: top corner left (column)
%   Column 2: top corner left (row)
%   Column 3: bottom corner right (column)
%   Column 4: bottom corner right (row)
b_mat = [Inf*ones(size(tracking_seq,1),1), Inf*ones(size(tracking_seq,1),1), ...
    -Inf*ones(size(tracking_seq,1),1), -Inf*ones(size(tracking_seq,1),1)];

for row_num = 1:size(tracking_seq,1)
    if (all(tracking_seq(row_num,:) <= 0))
        continue;
    end
    
    all_ad_overlap = zeros(size(image_sets{1}.objects));
    
    for col_num = 1:size(tracking_seq,2)
        if (tracking_seq(row_num,col_num) <= 0)
            continue;
        end
        
        this_puncta = image_sets{col_num}.objects == tracking_seq(row_num,col_num);
        all_ad_overlap = all_ad_overlap | this_puncta;
    end
    
    corners = [find(sum(all_ad_overlap),1,'first'), find(sum(all_ad_overlap,2),1,'first') ...
        find(sum(all_ad_overlap),1,'last'), find(sum(all_ad_overlap,2),1,'last')];
    b_mat(row_num,:) = corners;
end

image_size = size(image_sets{1}.objects);
b_mat(:,1:2) = b_mat(:,1:2) - i_p.Results.image_padding;
b_mat(:,1:2) = floor(b_mat(:,1:2));
b_mat(b_mat(:,1) <= 0,1) = 1;
b_mat(b_mat(:,2) <= 0,2) = 1;

b_mat(:,3:4) = b_mat(:,3:4) + i_p.Results.image_padding;
b_mat(:,3:4) = ceil(b_mat(:,3:4));
b_mat(b_mat(:,3) > image_size(2),3) = image_size(2);
b_mat(b_mat(:,4) > image_size(1),4) = image_size(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Single Ad Image Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puncta_seq = cell(size(tracking_seq,1), 1);
gel_seq = cell(size(tracking_seq,1), 1);

puncta_seq_high = cell(size(tracking_seq,1), 1);
gel_seq_high = cell(size(tracking_seq,1), 1);

for row_num = 1:size(tracking_seq,1)
    %some rows will be have been completely blanked by the above filtering,
    %skip these rows
    if (all(tracking_seq(row_num,:) <= 0))
        continue;
    end
    
    %we know we have a row with an object now, scan through each image,
    %pulling out that object
    for col_num = 1:size(tracking_seq,2)
        %attempt to find the surrounding entries of the current tracking
        %number position, if any of these numbers are above zero,
        %indicating that an object is present in this image or in the next
        %or immediately prior image
        surrounding_entries = [0, tracking_seq(row_num,col_num), 0];
        try surrounding_entries(1) = tracking_seq(row_num,col_num-1); end %#ok<TRYNC>
        try surrounding_entries(3) = tracking_seq(row_num,col_num+1); end %#ok<TRYNC>
        
        if (not(any(surrounding_entries > 0)))
            continue;
        end
        
        puncta_num = tracking_seq(row_num,col_num);
        %puncta_num of 0 indicates no object, but messes up the fact that
        %the background of the label matrices is 0, so reset the number to
        %NaN, so that their aren't any matches for highlighting before the
        %puncta of interest is present in the frame
        if (puncta_num == 0)
            puncta_num = NaN;
        end
        
        this_puncta_images = struct;
        this_puncta_images.puncta = image_sets{col_num}.objects == puncta_num;
        this_puncta_images.puncta_perim = image_sets{col_num}.objects_perim == puncta_num;
        this_puncta_images.not_this_puncta = image_sets{col_num}.objects ~= puncta_num & ...
            image_sets{col_num}.objects ~= 0;
        this_puncta_images.not_this_puncta_perim = image_sets{col_num}.objects_perim ~= puncta_num & ...
            image_sets{col_num}.objects_perim ~= 0;
        this_puncta_images.gel_image_norm = image_sets{col_num}.gel_image_norm;
        this_puncta_images.puncta_image_norm = image_sets{col_num}.puncta_image_norm;
        
        %remove all the area outside the bounds plus the padding of the
        %objects range based on the above bounding calculations
        f_names = fieldnames(this_puncta_images);
        for i = 1:length(f_names)
            this_field = f_names{i};
            this_puncta_images.(this_field) = bound_image(this_puncta_images.(this_field),b_mat(row_num,:));
        end
        
        background_region = logical(imdilate(this_puncta_images.puncta,strel('disk',5,0)));
        background_region = background_region & ...
            not(this_puncta_images.not_this_puncta) & ...
            not(this_puncta_images.puncta);
        
        %save a copy of the un-highlighted puncta and gel images for output
        %later
        puncta_seq{row_num}{col_num} = this_puncta_images.puncta_image_norm;
        gel_seq{row_num}{col_num} = this_puncta_images.gel_image_norm;
        
        %apply the highlights to the image sets
        highlight_fields = {'gel_image_norm','puncta_image_norm'};
        for i = 1:length(highlight_fields)
            this_field = highlight_fields{i};
            this_puncta_images.(this_field) = create_highlighted_image(this_puncta_images.(this_field), ...
                this_puncta_images.puncta_perim,'color_map',[0,1,0]);
            this_puncta_images.(this_field) = create_highlighted_image(this_puncta_images.(this_field), ...
                this_puncta_images.not_this_puncta_perim,'color_map',[0,0,1]);
            this_puncta_images.(this_field) = create_highlighted_image(this_puncta_images.(this_field), ...
                background_region,'color_map',[186/255,94/255,187/255],'mix_percent',0.5);
        end
        puncta_seq_high{row_num}{col_num} = this_puncta_images.puncta_image_norm;
        gel_seq_high{row_num}{col_num} = this_puncta_images.gel_image_norm;
    end
    
    puncta_seq{row_num} = remove_empty_cells(puncta_seq{row_num});
    gel_seq{row_num} = remove_empty_cells(gel_seq{row_num});
    puncta_seq_high{row_num} = remove_empty_cells(puncta_seq_high{row_num});
    gel_seq_high{row_num} = remove_empty_cells(gel_seq_high{row_num});
    
    if (i_p.Results.sample_number ~= 0)
        sample_indexes = round(linspace(1,length(puncta_seq{row_num}),i_p.Results.sample_number));
        puncta_seq_high{row_num} = puncta_seq_high{row_num}(sample_indexes);
        gel_seq_high{row_num} = gel_seq_high{row_num}(sample_indexes);
    end
    
    %Montage creation
    puncta_montage = create_montage_image_set(puncta_seq{row_num},'num_rows',1);
    gel_montage = create_montage_image_set(gel_seq{row_num},'num_rows',1);
    
    puncta_montage_high = create_montage_image_set(puncta_seq_high{row_num},'num_rows',1);
    gel_montage_high = create_montage_image_set(gel_seq_high{row_num},'num_rows',1);
    
    spacer = 0.5*ones(1,size(puncta_montage,2),3);
    spacer_high = 0.5*ones(1,size(puncta_montage_high,2),3);
    
    %Image output
    if (exist('invado_nums','var') && any(row_num == invado_nums))
        sub_out_folder = 'single_invadopodia';
    end
    if (exist('not_invado_nums','var') && any(row_num == not_invado_nums))
        sub_out_folder = 'single_notinvadopodia';
    end
    
    output_folder_no_high = fullfile(exp_dir,'visualizations',sub_out_folder,'no_highlights');
    if (not(exist(output_folder_no_high,'dir')))
        mkdir(output_folder_no_high);
    end
    
    output_folder_high = fullfile(exp_dir,'visualizations',sub_out_folder,'highlights');
    if (not(exist(output_folder_high,'dir')))
        mkdir(output_folder_high);
    end
    
    no_highlights = [puncta_montage;spacer;gel_montage];
    highlights = [puncta_montage_high;spacer_high;gel_montage_high];
    if (not(any(strcmp(i_p.UsingDefaults,'pixel_size'))))
        no_highlights = draw_scale_bar(no_highlights,i_p.Results.pixel_size,'bar_size',5);
        highlights = draw_scale_bar(highlights,i_p.Results.pixel_size,'bar_size',5);
    end
    
    imwrite(no_highlights, ...
        fullfile(output_folder_no_high,[num2str(row_num),'.png']));
    imwrite(highlights, ...
        fullfile(output_folder_high,[num2str(row_num),'.png']));
end
toc(start_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bounded_image = bound_image(image,bounding_mat)

bounded_image = image(bounding_mat(2):bounding_mat(4),bounding_mat(1):bounding_mat(3),:);

function image_set = remove_empty_cells(image_set)

while isempty(image_set{1}), image_set = image_set(2:end); end
while isempty(image_set{end}), image_set = image_set(1:(end-1)); end

function make_single_puncta_frames(exp_dir,varargin)
%MAKE_SINGLE_AD_FRAMES    Builds single image montages that track single
%                         adhesions through their entire lifecycle,
%                         including frames immediately preceding and
%                         following the adhesion's lifetime, if available
%
%   make_single_ad_frames(cfg_file,options) builds single adhesion montages
%   from raw experimental data, most config options are set in cfg_file
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'MAKE_SINGLE_AD_FRAMES';

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('start_row',1,@(x)x >= 1);
i_p.addParamValue('end_row',NaN,@(x)x >= 1);
i_p.addParamValue('puncta_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('image_padding',15,@isnumeric);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('..'));
filenames = add_filenames_to_struct(struct());
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking Sequence Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tracking_seq = load(fullfile(base_dir, image_dirs(1).name,filenames.tracking_matrix)) + 1;
% tracking_seq = tracking_seq(:,1:20);
% image_dirs = image_dirs(1:20);

start_row = i_p.Results.start_row;
end_row = i_p.Results.end_row;
if (isnan(end_row) || end_row > size(tracking_seq,1))
    end_row = size(tracking_seq,1);
end
tracking_temp = zeros(size(tracking_seq,1),size(tracking_seq,2));
tracking_temp(start_row:end_row,:) = tracking_seq(start_row:end_row,:);
tracking_seq = tracking_temp;

if (not(any(strcmp('puncta_file',i_p.UsingDefaults))))
    %the first row of the puncta data file is a header line
    puncta_data = csvread(i_p.Results.puncta_file,1,0);
    if (all(size(puncta_data) == [0 1]))
        puncta_data = 0;
    end
    
    puncta_nums = puncta_data(2:end,1);
    tracking_temp = zeros(size(tracking_seq,1),size(tracking_seq,2));
    tracking_temp(puncta_nums,:) = tracking_seq(puncta_nums,:);
    tracking_seq = tracking_temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Bounding Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The bounding matrix will hold the proper ranges to cut out only the
%relavent area for each adhesion's life cycle. The matrix is initialized
%with values that will always be replaced and the data is structured like
%this:
%   Column 1: top corner left (column)
%   Column 2: top corner left (row)
%   Column 3: bottom corner right (column)
%   Column 4: bottom corner right (row)
b_mat = [Inf*ones(size(tracking_seq,1),1), Inf*ones(size(tracking_seq,1),1), ...
    -Inf*ones(size(tracking_seq,1),1), -Inf*ones(size(tracking_seq,1),1)];

image_sets = cell(size(image_dirs,1),1);
for i = 1:size(image_dirs,1)
    image_sets{i} = read_in_file_set(fullfile(base_dir,image_dirs(i).name),filenames);
    if (mod(i,10) == 0)
        disp(['Finished Reading ', num2str(i), '/',num2str(size(image_dirs,1))]);
    end
end


for row_num = 1:size(tracking_seq,1)
    if (all(tracking_seq(row_num,:) <= 0))
        continue;
    end
    
    all_ad_overlap = zeros(size(image_sets{1}.puncta));
    
    for col_num = 1:size(tracking_seq,2)
        if (tracking_seq(row_num,col_num) <= 0)
            continue;
        end
        
        this_puncta = image_sets{col_num}.puncta == tracking_seq(row_num,col_num);
        all_ad_overlap = all_ad_overlap | this_puncta;
    end
    
    bounds = regionprops(all_ad_overlap,'BoundingBox');
    corners = [bounds(1).BoundingBox(1), bounds(1).BoundingBox(2)];
    corners = [corners, corners + bounds(1).BoundingBox(3:4)]; %#ok<AGROW>
    b_mat(row_num,:) = corners;
end

image_size = size(image_sets{1}.puncta);
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
    if (all(tracking_seq(row_num,:) <= 0))
        continue;
    end
    
    for col_num = 1:size(tracking_seq,2)
        surrounding_entries = [0, tracking_seq(row_num,col_num), 0];
        try surrounding_entries(1) = tracking_seq(row_num,col_num-1); end %#ok<TRYNC>
        try surrounding_entries(3) = tracking_seq(row_num,col_num+1); end %#ok<TRYNC>
        
        if (not(any(surrounding_entries > 0)))
            continue;
        end
        
        puncta_num = tracking_seq(row_num,col_num);
        if (puncta_num == 0)
            puncta_num = NaN;
        end
        
        this_set = struct;
        
        this_set.puncta = image_sets{col_num}.puncta_perim == puncta_num;
        this_set.not_this_puncta = image_sets{col_num}.puncta_perim ~= puncta_num & ...
            image_sets{col_num}.puncta_perim ~= 0;
        this_set.gel_image_norm = image_sets{col_num}.gel_image_norm;
        this_set.puncta_image_norm = image_sets{col_num}.puncta_image_norm;
        
        f_names = fieldnames(this_set);
        for i = 1:length(f_names)
            this_field = f_names{i};
            this_set.(this_field) = bound_image(this_set.(this_field),b_mat(row_num,:));
        end
        puncta_seq{row_num}{col_num} = this_set.puncta_image_norm;
        gel_seq{row_num}{col_num} = this_set.gel_image_norm;
        
        highlight_fields = {'gel_image_norm','puncta_image_norm'};
        for i = 1:length(highlight_fields)
            this_field = highlight_fields{i};
            this_set.(this_field) = create_highlighted_image(this_set.(this_field), ...
                this_set.puncta,'color_map',[0,1,0],'mix_percent',0.5);
            this_set.(this_field) = create_highlighted_image(this_set.(this_field), ...
                this_set.not_this_puncta,'color_map',[0,0,1],'mix_percent',0.5);
        end
        puncta_seq_high{row_num}{col_num} = this_set.puncta_image_norm;
        gel_seq_high{row_num}{col_num} = this_set.gel_image_norm;
    end
    
    puncta_montage = create_montage_image_set(puncta_seq{row_num},'num_rows',1);
    gel_montage = create_montage_image_set(gel_seq{row_num},'num_rows',1);
    
    puncta_montage_high = create_montage_image_set(puncta_seq_high{row_num},'num_rows',1);
    gel_montage_high = create_montage_image_set(gel_seq_high{row_num},'num_rows',1);
    
    spacer = 0.5*ones(1,size(puncta_montage,2),3);
    
    sub_out_folder = 'single_invadopodia';
    if (regexp(i_p.Results.puncta_file,'not_invado_data.csv'))
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
    
    imwrite([puncta_montage;spacer;gel_montage], ... 
        fullfile(output_folder_no_high,[num2str(row_num),'.png']));
    imwrite([puncta_montage_high;spacer;gel_montage_high], ... 
        fullfile(output_folder_high,[num2str(row_num),'.png']));
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bounded_image = bound_image(image,bounding_mat)

bounded_image = image(bounding_mat(2):bounding_mat(4),bounding_mat(1):bounding_mat(3),:);
function find_future_degrade_values(base_dir,tracking_mat,lineage_and_start,output_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'FIND_FUTURE_DEGRAD_VALUES';

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('tracking_mat',@(x)exist(x,'file') == 2);
i_p.addRequired('lineage_and_start',@(x)exist(x,'file') == 2);
i_p.addRequired('output_file',@ischar);

i_p.addParamValue('adhesions_filename','puncta_labeled.png',@ischar);
i_p.addParamValue('gel_filename','registered_gel.png',@ischar);
i_p.addParamValue('binary_shift_filename','binary_shift.png',@ischar);

% i_p.addParamValue('output_dir', fileparts(adhesions_file), @(x)exist(x,'dir')==7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(base_dir,tracking_mat,lineage_and_start,output_file,varargin{:});

tracking_mat = csvread(i_p.Results.tracking_mat);
lineage_and_start = csvread(i_p.Results.lineage_and_start);
lineage_nums = lineage_and_start(1,:);
start_nums = lineage_and_start(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% base_dir = fileparts(fileparts(i_p.Results.adhesions_file));
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

% image_path = fileparts(i_p.Results.adhesions_file);
% [this_path, this_image_num, ext, version] = fileparts(image_path); %#ok<NASGU>

Local_gel_mat = NaN*ones(size(tracking_mat));

for i = 1:length(lineage_nums)
    adhesions = imread(fullfile(base_dir, image_dirs(start_nums(i)).name, i_p.Results.adhesions_filename));
    
    ad_num = tracking_mat(lineage_nums(i), start_nums(i)) + 1;
    adhesions_temp = zeros(size(adhesions));
    adhesions_temp(adhesions == ad_num) = 1;
    adhesions = adhesions_temp;
    
    post_death_cols = find(tracking_mat(lineage_nums(i),:) <= -1);
    post_death_cols(post_death_cols <= start_nums(i)) = 0;
    post_death_cols = unique(post_death_cols);
    assert(post_death_cols(1) == 0);
    assert(post_death_cols(2) > 0);
    post_death_cols = post_death_cols(2:end);
    
    for j = post_death_cols
        current_adhesions = imread(fullfile(base_dir, image_dirs(j).name, i_p.Results.adhesions_filename));
        binary_shift = logical(imread(fullfile(base_dir, image_dirs(j).name, i_p.Results.binary_shift_filename)));
        gel = imread(fullfile(base_dir, image_dirs(j).name, i_p.Results.gel_filename));
        scale_factor = double(intmax(class(gel)));
        gel = double(gel)/scale_factor;
        
        temp = find_local_gel_diff(adhesions,gel,binary_shift, current_adhesions);
        
        Local_gel_mat(lineage_nums(i),j) = temp(1).Local_gel_diff;
    end
    disp(['Done with ', num2str(lineage_nums(i))]);
end
disp('Done with determining local diff.');

% dlmwrite(fullfile(i_p.Results.output_dir,'raw_data','Local_gel_mat.csv'), Local_gel_mat,'precision','%.4f')
dlmwrite(output_file, Local_gel_mat,'precision','%.8f')
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = find_local_gel_diff(labeled_adhesions,gel_image,binary_shift,current_adhesions,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_LOCAL_GEL_DIFF';

i_p.addRequired('gel_image',@isnumeric);
i_p.addRequired('labeled_adhesions',@(x)isnumeric(x));
i_p.addRequired('binary_shift',@(x)islogical(x));
i_p.addRequired('current_adhesions',@(x)isnumeric(x));

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addParamValue('cell_mask',@(x)islogical(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(labeled_adhesions,gel_image,binary_shift,current_adhesions,varargin{:});

adhesion_props = regionprops(labeled_adhesions,'Area');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(labeled_adhesions(:))
    if (not(any(any(labeled_adhesions == i))))
        adhesion_props(i).Local_gel_diff = 0;
        adhesion_props(i).Global_gel_diff = 0;
        continue;
    end
    
    this_ad = labeled_adhesions;
    this_ad(labeled_adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
    background_region = and(background_region,not(labeled_adhesions));
    background_region = and(background_region,not(current_adhesions));
    background_region = and(background_region, binary_shift);
    
    assert(sum(sum(background_region)) > 0)
    
    adhesion_props(i).Local_gel_diff = mean(gel_image(this_ad)) - mean(gel_image(background_region));
    adhesion_props(i).Global_gel_diff = mean(gel_image(this_ad)) - mean(gel_image(binary_shift));
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(labeled_adhesions(:)))]); end
end
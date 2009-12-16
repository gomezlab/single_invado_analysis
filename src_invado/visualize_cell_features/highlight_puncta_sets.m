function highlight_puncta_sets(cfg_file,i_num,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'HIGHLIGHT_PUNCTA_SETS';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addRequired('i_num',@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,i_num,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(cfg_file);
while 1
    line = fgetl(fid);
    if ~ischar(line), break; end
    eval(line);
end

addpath(genpath(path_folders));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect General Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_image_num = find_max_image_num(I_folder);
folder_char_length = length(num2str(max_image_num));
i_size = size(imread(fullfile(I_folder,num2str(max_image_num),focal_image)));

tracking_seq = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find all the image directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_directories = dir(I_folder);

assert(not(isempty(strmatch(image_directories(1).name,'.','exact'))), ... 
    'Expected first entry in image_directories to be "."');
assert(not(isempty(strmatch(image_directories(2).name,'..','exact'))), ... 
    'Expected second entry in image_directories to be ".."');

image_directories = image_directories(3:end);

for i=1:size(image_directories,1)
    assert(str2num(image_directories(i).name) == i) %#ok<ST2NM>
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in all the files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puncta_perim = imread(fullfile(I_folder,image_directories(i_num).name,adhesions_perim_filename));

binary_shift = imread(fullfile(I_folder,image_directories(i_num).name,'binary_shift.png'));

puncta_image = imread(fullfile(I_folder,image_directories(i_num).name,focal_image));
puncta_image_range = csvread(fullfile(I_folder,image_directories(i_num).name,'puncta_image_range.csv'));
puncta_image = normalize_image(puncta_image,binary_shift,puncta_image_range);
orig_puncta = puncta_image;

gel_image = imread(fullfile(I_folder,image_directories(i_num).name,'registered_gel.png'));
gel_image_range = csvread(fullfile(I_folder,image_directories(i_num).name,'gel_image_range.csv'));
gel_image = normalize_image(gel_image,binary_shift,gel_image_range);
orig_gel = gel_image;

invado_nums = csvread(fullfile(lin_time_series_folder, '../invado_nums.csv'));
non_invado_nums = csvread(fullfile(lin_time_series_folder, '../non_invado_nums.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Highlight the appropriate adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puncta_image = create_highlighted_image(puncta_image,im2bw(puncta_perim,0),'color_map',[0,0,1]);
gel_image = create_highlighted_image(gel_image,im2bw(puncta_perim,0),'color_map',[0,0,1]);

invado_tracking = tracking_seq(invado_nums,:);
non_invado_tracking = tracking_seq(non_invado_nums,:);

invado_image_nums = sort(unique(invado_tracking(:,i_num)));
non_invado_image_nums = sort(unique(non_invado_tracking(:,i_num)));

if (not(all(invado_image_nums <= 0)))
    invado_image_nums = remove_neg_nums(invado_image_nums);
    
    highlight_binary = zeros(i_size);
    for i = 1:size(invado_image_nums,1)
        highlight_binary(puncta_perim == invado_image_nums(i)) = 1;
    end
    
    puncta_image = create_highlighted_image(puncta_image,highlight_binary,'color_map',[0,1,0]);
    gel_image = create_highlighted_image(gel_image,highlight_binary,'color_map',[0,1,0]);
end

if (not(all(non_invado_image_nums <= 0)))
    non_invado_image_nums = remove_neg_nums(non_invado_image_nums);
    
    highlight_binary = zeros(i_size);
    for i = 1:size(non_invado_image_nums,1)
        highlight_binary(puncta_perim == non_invado_image_nums(i)) = 1;
    end
    
    puncta_image = create_highlighted_image(puncta_image,highlight_binary,'color_map',[1,0,0]);
    gel_image = create_highlighted_image(gel_image,highlight_binary,'color_map',[1,0,0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Composite and Output the Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dual_image = [gel_image, 0.5*ones(size(puncta_image,1),1,3), puncta_image];

imwrite(dual_image, fullfile(I_folder,image_directories(i_num).name,'invado_and_not.png'));

function normed_image = normalize_image(image, binary_shift, range)

normed_image = (double(image) - range(1))*(1/range(2));
normed_image(not(binary_shift)) = 0;


function no_negatives = remove_neg_nums(this_matrix)

assert(size(this_matrix,1) == 1 || size(this_matrix,2) == 1);
sorted_mat = sort(this_matrix);
first_pos_index = find(sorted_mat > 0,1,'first');
no_negatives = sorted_mat(first_pos_index:end);
assert(all(no_negatives > 0))
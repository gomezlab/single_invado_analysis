function find_future_degrade_values(adhesions_file,gel_file,binary_shift_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'FIND_FUTURE_DEGRAD_VALUES';

i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);
i_p.addRequired('gel_file',@(x)exist(x,'file') == 2);
i_p.addRequired('binary_shift_file',@(x)exist(x,'file') == 2);

i_p.addParamValue('output_dir', fileparts(adhesions_file), @(x)exist(x,'dir')==7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(adhesions_file,gel_file,binary_shift_file,varargin{:});

%determine the file name used for the gel files
gel_filename = dir(gel_file);

%read in the labeled adhesions
adhesions = imread(adhesions_file);
adhesions_filename = dir(adhesions_file);

%determine the file name used for the binary shift files
binary_shift_filename = dir(binary_shift_file);

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fileparts(fileparts(i_p.Results.adhesions_file));
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

image_path = fileparts(i_p.Results.adhesions_file);
[this_path, this_image_num, ext, version] = fileparts(image_path); %#ok<NASGU>

Local_gel_mat = zeros(max(adhesions(:)), length(image_dirs) - (str2num(this_image_num)+1)); %#ok<ST2NM>
gel_increment = 0;

for i=(str2num(this_image_num)+1):length(image_dirs) %#ok<ST2NM>
    current_adhesions = imread(fullfile(base_dir, image_dirs(i).name, adhesions_filename.name));
    binary_shift = logical(imread(fullfile(base_dir, image_dirs(i).name, binary_shift_filename.name)));
    gel = imread(fullfile(base_dir, image_dirs(i).name, gel_filename.name));
    
    temp = find_local_gel_diff(adhesions,gel,binary_shift, current_adhesions);
    gel_increment = gel_increment + 1;
    Local_gel_mat(:,gel_increment) = [temp.Local_gel_diff];
    
    disp(['Finished: ',fullfile(base_dir, image_dirs(i).name)])
end

dlmwrite(fullfile(i_p.Results.output_dir,'raw_data','Local_gel_mat.csv'), Local_gel_mat,'precision','%.4f')

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
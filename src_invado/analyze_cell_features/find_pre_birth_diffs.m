function find_pre_birth_diffs(cfg_file,image_num,varargin)
%FIND_PRE_BIRTH_DIFFS    Searches through a given tracking matrix and a
%                        data set to produce local diff values for each
%                        puncta immediately before the puncta's birth, if
%                        such a time is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_PRE_BIRTH_DIFFS';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addRequired('image_num',@(x)isnumeric(x) && x > 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.addParamValue('adhesions_filename','puncta_labeled.png',@ischar);
i_p.addParamValue('puncta_filename','registered_focal_image.png',@ischar);
i_p.addParamValue('gel_filename','registered_gel.png',@ischar);
i_p.addParamValue('binary_shift_filename','binary_shift.png',@ischar);
i_p.addParamValue('cell_mask_filename','cell_mask.png',@ischar);
i_p.addParamValue('intensity_correction_file','intensity_correction.csv',@ischar);

i_p.parse(cfg_file,image_num,varargin{:});

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
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%after loading the tracking sequence filter to only include those puncta
%included in the invadopodia list, remember the list formated so the first
%column contains the lineage number, with the first lineage as 1, so no
%need to translate
tracking_seq = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%final all the image directories
image_dirs = dir(I_folder);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);
assert(image_num <= length(image_dirs));

ad_label = imread(fullfile(I_folder,image_dirs(image_num).name,adhesions_filename));

pre_birth_diffs = ones(1,max(ad_label(:)))*NaN;
gel_limits = csvread(fullfile(I_folder,image_dirs(end).name,'gel_image_range.csv'));

for i = 1:size(tracking_seq,1)
    %there won't ever be pre birth data when the image number is 1, as the
    %puncta detected in image 1 don't have a birth event
    if (image_num == 1)
        continue;
    end
    
    %skip over any rows that don't have a real puncta number
    if (tracking_seq(i,image_num) <= 0)
        continue;
    end
    
    %now we have a real puncta number, lets find the pre-birth image number
    birth_i_num = find_birth_i_num(tracking_seq(i,:) < 1) - 1;
    %the above function returns NaN if there isn't a pre-birth image
    %number, catch that situation and exit out of the loop in that case
    if (isnan(birth_i_num))
        continue;
    end
    
    image_data = struct();
    
    image_data.binary_shift = logical(imread(fullfile(I_folder,image_dirs(birth_i_num).name,i_p.Results.binary_shift_filename)));
    
    %read in the gel image and normalize to 0-1
    image_data.gel_image = double(imread(fullfile(I_folder,image_dirs(birth_i_num).name,i_p.Results.gel_filename)));
    image_data.gel_image = image_data.gel_image - gel_limits(1);
    image_data.gel_image = image_data.gel_image .* (1/gel_limits(2));
    image_data.gel_image(not(image_data.binary_shift)) = 0;
    image_data.gel_image = cat(3,image_data.gel_image,image_data.gel_image,image_data.gel_image);
    
    image_data.intensity_correction = csvread(fullfile(I_folder,image_dirs(birth_i_num).name,i_p.Results.intensity_correction_file));
    image_data.adhesions = imread(fullfile(I_folder,image_dirs(birth_i_num).name, i_p.Results.adhesions_filename));
    
    puncta_num = tracking_seq(i,image_num);
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = ad_label;
    this_ad(ad_label ~= puncta_num) = 0;
    this_ad = logical(this_ad);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    this_pre_birth_diff = collect_local_diff_properties(image_data,this_ad);
    
    pre_birth_diffs(puncta_num) = this_pre_birth_diff.Local_gel_diff;
    
    if(i_p.Results.debug), disp([puncta_num,birth_i_num]); end
end

dlmwrite(fullfile(I_folder,image_dirs(image_num).name,'raw_data','Pre_birth_diff.csv'), pre_birth_diffs');

profile off;
if (i_p.Results.debug), profile viewer; end
end

function birth_i_num = find_birth_i_num(puncta_present_logical)
%FIND_BIRTH_I_NUM    Searches through a logical matrix for switches from
%                    positive to negative and vice versa that indicate
%                    birth and death events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

birth_i_num = NaN;
transistion_points = [];
for i=2:length(puncta_present_logical)
    if (xor(puncta_present_logical(i), puncta_present_logical(i-1)))
        transistion_points = [transistion_points, i]; %#ok<AGROW>
    end
end

assert(length(transistion_points) <= 2);

if (length(transistion_points) >= 1)
    birth_i_num = transistion_points(1);
end
end
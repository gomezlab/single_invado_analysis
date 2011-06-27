function find_pre_birth_diffs(exp_folder,image_num,varargin)
%FIND_PRE_BIRTH_DIFFS    Searches through a given tracking matrix and a
%                        data set to produce local diff values for each
%                        puncta immediately before the puncta's birth, if
%                        such a time is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'FIND_PRE_BIRTH_DIFFS';

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);
i_p.addRequired('image_num',@(x)isnumeric(x) && x > 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_folder,image_num,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = fullfile(exp_folder,'individual_pictures');

%find all the image directories
image_dirs = dir(base_folder);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);
assert(image_num <= length(image_dirs));

tracking_seq = load(fullfile(base_folder,image_dirs(1).name,filenames.tracking_matrix))+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(tracking_seq,2)
    data_set{j} = read_in_file_set(fullfile(base_folder,image_dirs(j).name),filenames); %#ok<AGROW>
    pre_birth_diffs{j} = ones(1,max(data_set{j}.puncta(:)))*NaN; %#ok<NASGU>
end

current_data = data_set{image_num};

pre_birth_diffs = ones(1,max(current_data.puncta(:)))*NaN;

for row = 1:size(tracking_seq,1)
    %there won't ever be pre birth data when the image number is 1, as the
    %puncta detected in image 1 don't have a birth event
    if (image_num == 1)
        continue;
    end
    
    %skip over any rows that don't have any puncta
    if (tracking_seq(row,image_num) <= 0)
        continue;
    end
    
    %now we have a real puncta number, lets find the pre-birth image number
    pre_birth_i_num = find_birth_i_num(tracking_seq(row,:) < 1) - 1;
    %the above function returns NaN if there isn't a pre-birth image
    %number, catch that situation and exit out of the loop in that case
    if (isnan(pre_birth_i_num))
        continue;
    end
    
    pre_birth_image_data = data_set{j};
    
    read_in_file_set(fullfile(base_folder,image_dirs(pre_birth_i_num).name),filenames);
        
    puncta_num = tracking_seq(row,image_num);
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = current_data.puncta;
    this_ad(current_data.puncta ~= puncta_num) = 0;
    this_ad = logical(this_ad);
        
    this_pre_birth_diff = collect_local_diff_properties(pre_birth_image_data,this_ad);
    
    pre_birth_diffs(puncta_num) = this_pre_birth_diff.Local_gel_diff;
    
    if(i_p.Results.debug), disp([puncta_num,pre_birth_i_num]); end
end

dlmwrite(fullfile(base_folder,image_dirs(image_num).name,'raw_data','Pre_birth_diff.csv'), pre_birth_diffs');

toc;
end

function birth_i_num = find_birth_i_num(puncta_present_logical)
%FIND_BIRTH_I_NUM    Searches through a logical matrix for switches from
%                    positive to negative and vice versa that indicate
%                    birth and death events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

birth_i_num = NaN;
for i=2:length(puncta_present_logical)
    if (puncta_present_logical(i) == 0 && puncta_present_logical(i-1) == 1)
        %if the birth_i_num is not nan, we have an error since we should
        %only see one transition from 1 to 0
        assert(isnan(birth_i_num))
            
        birth_i_num = i;
    end
end

end

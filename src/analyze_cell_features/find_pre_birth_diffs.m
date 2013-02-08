function find_pre_birth_diffs(exp_folder,varargin)
%FIND_PRE_BIRTH_DIFFS    Searches through a given tracking matrix and a
%                        data set to produce local diff values for each
%                        puncta immediately before the puncta's birth, if
%                        such a time is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_all=tic;
i_p = inputParser;
i_p.FunctionName = 'FIND_PRE_BIRTH_DIFFS';

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_folder,varargin{:});

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

tracking_seq = load(fullfile(base_folder,image_dirs(1).name,filenames.tracking_matrix))+1;
pre_birth_diffs = NaN*ones(size(tracking_seq));

% tracking_seq = tracking_seq(:,1:20);
% image_dirs = image_dirs(1:20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read In Data Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_read = tic;
for j = 1:size(tracking_seq,2)
    data_set{j} = read_in_file_set(fullfile(base_folder,image_dirs(j).name),filenames); %#ok<AGROW>
    if (mod(j,10) == 0)
        disp(['Read in ',num2str(j),'/',num2str(size(tracking_seq,2))]);
    end
end
toc(start_read);

start_process = tic;
for row = 1:size(tracking_seq,1)
    pre_birth_i_num = find_birth_i_num(tracking_seq(row,:) < 1) - 1;
    %the above function returns NaN if there isn't a pre-birth image
    %number, catch that situation and exit out of the loop in that case
    if (isnan(pre_birth_i_num))
        continue;
    end
    
    pre_birth_image_data = data_set{pre_birth_i_num};
    
    for col = 1:size(tracking_seq,2)
        if (tracking_seq(row,col) <= 0)
            continue;
        end
        
        puncta_num = tracking_seq(row,col);
        
        this_puncta = data_set{col}.objects == puncta_num;
        pre_birth_diff = collect_local_diff_properties(pre_birth_image_data,this_puncta);
        
        pre_birth_diffs(row,col) = pre_birth_diff.Local_gel_diff;
    end
    
    if (mod(row,100) == 0)
        disp(['Processed ',num2str(row),'/',num2str(size(tracking_seq,1))]);
    end    
end
toc(start_process);

output_folder = fullfile(exp_folder,'puncta_props','lin_time_series');
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end
dlmwrite(fullfile(output_folder,'Pre_birth_diff.csv'), pre_birth_diffs);

toc(start_all);
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

function track_invado(exp_folder,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_all=tic;
i_p = inputParser;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Puncta Matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_match = tic;
matching_matrices = cell(size(image_dirs,1),1);
puncta_nums = cell(size(image_dirs,1),1);
tracking_counts = struct ('recip_best',0,'area',0);
for j = 1:size(image_dirs,1)
    puncta_1 = imread(fullfile(base_folder,image_dirs(j).name,filenames.objects));
    
    %the last image set is a special case, there won't be any more matches,
    %but we need place holders in the match sets, so that they get added to
    %final tracking matrix
    if (j == size(image_dirs,1))
        all_match_sets = zeros(length(puncta_nums{j}),2);
        matching_matrices{j} = all_match_sets;
        continue;
    end
        
    puncta_2 = imread(fullfile(base_folder,image_dirs(j+1).name,filenames.objects));
    puncta_nums{j} = 1:max(puncta_1(:));
    puncta_nums{j+1} = 1:max(puncta_2(:));

    counts = determine_pixel_count_overlaps(puncta_1,puncta_2);

    all_match_sets = zeros(length(puncta_nums{j}),2);
    while (not(exist('old_match_sets','var')) || ...  
           any(any(old_match_sets ~= all_match_sets)))
        old_match_sets = all_match_sets;
        [all_match_sets,decision_count] = find_dual_max_hits_sets(counts,all_match_sets);
        tracking_counts.recip_best = tracking_counts.recip_best + decision_count;
        
        [puncta_1,puncta_2] = remove_match_set_puncta(puncta_1,puncta_2,all_match_sets);
        counts = determine_pixel_count_overlaps(puncta_1,puncta_2);
        1;
    end
    if (any(any(counts > 0)))
        [all_match_sets,decision_count] = find_area_row_matches(counts,all_match_sets,puncta_1,puncta_2);
        tracking_counts.area = tracking_counts.area + decision_count;
        [puncta_1,puncta_2] = remove_match_set_puncta(puncta_1,puncta_2,all_match_sets);
        counts = determine_pixel_count_overlaps(puncta_1,puncta_2);

        [all_match_sets,decision_count] = find_area_col_matches(counts,all_match_sets,puncta_1,puncta_2);
        tracking_counts.area = tracking_counts.area + decision_count;
        [puncta_1,puncta_2] = remove_match_set_puncta(puncta_1,puncta_2,all_match_sets);
        counts = determine_pixel_count_overlaps(puncta_1,puncta_2);
    end

    matching_matrices{j} = all_match_sets;
    clear all_match_sets;
    clear old_match_sets;
    if (mod(j,1) == 0)
        disp(['Tracked ',num2str(j),'/',num2str(size(image_dirs,1))]);
    end
end
tracking_counts = process_tracking_counts(tracking_counts);
disp(['Tracking done, made ', tracking_counts.total, ' decisions,' ...
    tracking_counts.recip_percent,'% tracked with best reciprocal match, '...
    tracking_counts.area_percent, '% with largest area match.']);
toc(start_match);

% save('matching_mat.mat','matching_matrices');
% save('puncta_nums.mat','puncta_nums');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matching to Tracking Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('matching_mat.mat'); load('puncta_nums.mat');

%run through all the matching matices, building all the rows of the
%tracking matrix
track_rows = struct;
track_num = 1;
for i_num = 1:size(matching_matrices,1)
    for row = 1:size(matching_matrices{i_num},1);
        this_i_num = i_num;
        match_row = matching_matrices{this_i_num}(row,:);

        %this row of the matching matrix has already been cleared, skip to
        %next
        if (match_row(1) <= 0)
            if (match_row(1) == 0)
                track_rows(track_num).start = i_num;
                track_rows(track_num).seq = row;
                track_num = track_num + 1;
            end
            continue;
        end
        
        track_rows(track_num).start = i_num;
        track_rows(track_num).seq = match_row;
        
        matching_matrices{this_i_num}(row,:) = [-1,-1];
        this_i_num = this_i_num + 1;
        next_match_row = match_row(2);
        
        while (this_i_num <= size(matching_matrices,1) && next_match_row ~= 0)
            match_row = matching_matrices{this_i_num}(next_match_row,:);
            matching_matrices{this_i_num}(next_match_row,:) = [-1,-1];
            next_match_row = match_row(2);
            if (next_match_row == 0)
                continue;
            end
            
            assert(track_rows(track_num).seq(end) == match_row(1));
            track_rows(track_num).seq = [track_rows(track_num).seq,match_row(2)];
            
            this_i_num = this_i_num + 1;
        end
        track_num = track_num + 1;
    end
    assert(all(all(matching_matrices{i_num} == 0 | matching_matrices{i_num} == -1)))
end

%convert the tracking rows to a tracking matrix
total_i_nums = size(image_dirs,1);
tracking_matrix = zeros(size(track_rows,2),total_i_nums);

for track_num = 1:size(track_rows,2)
    tracking_row = [zeros(1,track_rows(track_num).start-1), ...
        track_rows(track_num).seq, ...
        zeros(1,length(image_dirs)-length(track_rows(track_num).seq)-(track_rows(track_num).start-1))];
    
    tracking_row = [tracking_row,zeros(1,total_i_nums - length(tracking_row))]; %#ok<AGROW>
    assert(length(tracking_row) == total_i_nums);
    
    tracking_matrix(track_num,:) = tracking_row;
end

%check each column of the tracking matrix to make sure every object in the
%original images are present in the tracking matrix
for col_num = 1:size(tracking_matrix,2)
    these_puncta_nums = unique(tracking_matrix(:,col_num));
    if (these_puncta_nums(1) == 0), these_puncta_nums = these_puncta_nums(2:end); end
    
    assert(all(these_puncta_nums == puncta_nums{col_num}'));
end

csvwrite(fullfile(exp_folder,'tracking_matrices','tracking_seq_diff.csv'),tracking_matrix - 1);
toc(start_all);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best Match Finding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function counts = determine_pixel_count_overlaps(puncta_1,puncta_2)
counts = zeros(max(puncta_1(:)),max(puncta_2(:)));
for puncta_num = 1:max(puncta_1(:))
    this_puncta = puncta_1 == puncta_num;
    
    overlap_puncta = puncta_2 .* uint16(this_puncta);
    
    overlap_nums = unique(overlap_puncta);
    if (overlap_nums(1) == 0), overlap_nums = overlap_nums(2:end); end
    if (isempty(overlap_nums)), continue; end
    
    for this_overlap_num = overlap_nums';
        counts(puncta_num,this_overlap_num) = sum(sum(overlap_puncta == this_overlap_num));
    end
end

end

function [dual_match_sets,decision_count] = find_dual_max_hits_sets(counts,dual_match_sets)

decision_count = 0;
for row = 1:size(counts,1)
    %if there isn't any overlap, skip out of the loop, there won't be any
    %recriprocal hits
    if (all(counts(row,:) == 0))
        continue;
    end
    
    %max_row_matches holds the column numbers that match the highest value
    %in the row
    max_row_matches = find(counts(row,:) == max(counts(row,:)));
    
    %these commands execute the reciprocal search, namely, what row matches
    %the highest value in the column
    col = max_row_matches(1);
    max_col_matches = find(counts(:,col) == max(counts(:,col)));
    
    %if the max rows and columns are unique and they reciprocally match, we
    %have a unique two way hit, save the hit
    if (length(max_row_matches) == 1 && length(max_col_matches) == 1 && ...
            col == max_row_matches(1) && row == max_col_matches(1))
        dual_match_sets(row,:) = [row,col];
        decision_count = decision_count + 1;
    end
end

end

function [cleared_puncta_1,cleared_puncta_2] = remove_match_set_puncta(puncta_1,puncta_2,match_sets)

cleared_puncta_1 = puncta_1;
cleared_puncta_2 = puncta_2;

for row = 1:size(match_sets,1)
    if (match_sets(row,1) == 0)
        continue;
    end
    cleared_puncta_1(puncta_1 == match_sets(row,1)) = 0;
    cleared_puncta_2(puncta_2 == match_sets(row,2)) = 0;
end

end

function [dual_match_sets,decision_count] = find_area_row_matches(counts,dual_match_sets,puncta_1,puncta_2)

decision_count = 0;
for row = 1:size(counts,1)
    %if there isn't any overlap, skip out of the loop, there won't be any
    %row matches
    if (all(counts(row,:) == 0))
        continue;
    end
    
    %max_row_matches holds the column numbers that match the highest value
    %in the row
    max_row_matches = find(counts(row,:) == max(counts(row,:)));
    
    if (length(max_row_matches) < 2)
        continue;
    end
    
    next_props = regionprops(puncta_2,'Area');
    
    largest_area = -Inf;
    largest_puncta_num = 0;
    
    for next_puncta_num = max_row_matches
        if (next_props(next_puncta_num).Area > largest_area)
            largest_area = next_props(next_puncta_num).Area;
            largest_puncta_num = next_puncta_num;
        end
    end
    
    dual_match_sets(row,:) = [row,largest_puncta_num];
    decision_count = decision_count + 1;
end

end

function [dual_match_sets,decision_count] = find_area_col_matches(counts,dual_match_sets,puncta_1,puncta_2)

decision_count = 0;
for col = 1:size(counts,2)
    %if there isn't any overlap, skip out of the loop, there won't be any
    %row matches
    if (all(counts(:,col) == 0))
        continue;
    end
    
    %max_row_matches holds the column numbers that match the highest value
    %in the row
    max_col_matches = find(counts(:,col) == max(counts(:,col)));
    
    if (length(max_col_matches) < 2)
        continue;
    end
    
    cur_props = regionprops(puncta_1,'Area');
    
    largest_area = -Inf;
    largest_puncta_num = 0;
    
    for next_puncta_num = max_col_matches'
        if (cur_props(next_puncta_num).Area > largest_area)
            largest_area = cur_props(next_puncta_num).Area;
            largest_puncta_num = next_puncta_num;
        end
    end
    
    dual_match_sets(largest_puncta_num,:) = [largest_puncta_num,col];
    decision_count = decision_count + 1;
end

end

function tracking_counts = process_tracking_counts(tracking_counts)

tracking_counts.total = tracking_counts.recip_best + tracking_counts.area;
tracking_counts.recip_percent = 100*tracking_counts.recip_best/tracking_counts.total;
tracking_counts.area_percent = 100*tracking_counts.area/tracking_counts.total;

tracking_counts.recip_percent = sprintf('%0.2f',tracking_counts.recip_percent);
tracking_counts.area_percent = sprintf('%0.2f',tracking_counts.area_percent);

end
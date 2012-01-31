function process_results(base_dir,property_folder,degrade_frames,cell_num,base_file_name)
% base_dir = '../../results/Invado_count_20000/';
% property_folder = '2011_12_02/control/time_series_01/cell_props/lin_time_series/';
% degrade_frames = [7 8 9 10];
% cell_num = 2;
% base_file_name = '2011_12_02-co-time_series_01-002';

full_property_folder = fullfile(base_dir,property_folder);
results = cell(0);

if (not(exist(fullfile(full_property_folder,'Cell_gel_diff.csv'),'file')))
    base_file_name
    return;
end

results.centroid_x = csvread(fullfile(full_property_folder,'Centroid_x.csv'));
results.centroid_x = results.centroid_x(cell_num,:);

results.centroid_y = csvread(fullfile(full_property_folder,'Centroid_y.csv'));
results.centroid_y = results.centroid_y(cell_num,:);

results.diff_p_val = csvread(fullfile(full_property_folder,'Cell_gel_diff_p_val.csv'));
results.diff_p_val = results.diff_p_val(cell_num,:);

results.diff_mean = csvread(fullfile(full_property_folder,'Cell_gel_diff.csv'));
results.diff_mean = results.diff_mean(cell_num,:);

results.diff_mean_no_corr = csvread(fullfile(full_property_folder,'Cell_gel_diff_no_corr.csv'));
results.diff_mean_no_corr = results.diff_mean_no_corr(cell_num,:);

results.diff_percent = csvread(fullfile(full_property_folder,'Cell_gel_diff_percent.csv'));
results.diff_percent = results.diff_percent(cell_num,:);

results.diff_median = csvread(fullfile(full_property_folder,'Cell_gel_diff_median.csv'));
results.diff_median = results.diff_median(cell_num,:);

results.degrade_or_not = ones(1,length(results.diff_mean));
results.degrade_or_not(setdiff(1:length(results.diff_mean),degrade_frames)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing non calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we can't call a degradation event without the frame before it present, so
%we need to remove all the non-callable values up to the second data
%instance

callable = find(not(isnan(results.centroid_x)));
callable = callable(2:end);

fields = fieldnames(results);
for i=1:length(fields)
    results.(fields{i}) = results.(fields{i})(callable);
end

save(fullfile('combined_results/',[base_file_name,'.mat']),'results');
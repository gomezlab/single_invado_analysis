function find_invading_cells(field_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('field_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('median_filter',-1.5,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(field_dir,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data = struct();
files = struct();

data_series_folder = fullfile(field_dir,'cell_props','lin_time_series');

files.p_vals = fullfile(data_series_folder,'Cell_gel_diff_p_val.csv');
files.cell_diffs = fullfile(data_series_folder,'Cell_gel_diff.csv');
files.cell_diff_medians = fullfile(data_series_folder,'Cell_gel_diff_median.csv');
files.cell_diff_percents = fullfile(data_series_folder,'Cell_gel_diff_percent.csv');
files.tracking = fullfile(field_dir,'tracking_matrices','tracking_seq.csv');

these_types = fieldnames(files);
for j = 1:length(these_types)
    this_file = files.(these_types{j});
    
    %matlab doesn't like you to reference fields that haven't been
    %created, so create files that aren't present yet before loading
    %data in
    if(isempty(strmatch(these_types{j},fieldnames(raw_data))))
        raw_data.(these_types{j}) = [];
    end
    
    if (exist(this_file,'file'))
        raw_data.(these_types{j}) = load(this_file);
    else
        error('Invado:MissingFile',['Can''t find ',this_file])
    end
end

%check that all the raw data files are the same size
these_names = fieldnames(raw_data);
poss_name_combinations = combnk(1:length(these_names),2);
for j=1:size(poss_name_combinations,1)
    assert(all(size(raw_data.(these_names{poss_name_combinations(j,1)})) == ...
        size(raw_data.(these_names{poss_name_combinations(j,2)}))))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
processed_data = process_raw_data(raw_data);

output_dir = fullfile(field_dir,'cell_props');

csvwrite(fullfile(output_dir,'active_degrade.csv'),processed_data.active_degrade);

csvwrite(fullfile(output_dir,'active_degrade_10.csv'),processed_data.active_degrade_10);
csvwrite(fullfile(output_dir,'active_degrade_20.csv'),processed_data.active_degrade_20);
csvwrite(fullfile(output_dir,'active_degrade_30.csv'),processed_data.active_degrade_30);
csvwrite(fullfile(output_dir,'active_degrade_40.csv'),processed_data.active_degrade_40);
csvwrite(fullfile(output_dir,'active_degrade_50.csv'),processed_data.active_degrade_50);
csvwrite(fullfile(output_dir,'active_degrade_60.csv'),processed_data.active_degrade_60);
csvwrite(fullfile(output_dir,'active_degrade_70.csv'),processed_data.active_degrade_70);
csvwrite(fullfile(output_dir,'active_degrade_80.csv'),processed_data.active_degrade_80);
csvwrite(fullfile(output_dir,'active_degrade_90.csv'),processed_data.active_degrade_90);

csvwrite(fullfile(output_dir,'longevity.csv'),processed_data.longevities);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function process_data = process_raw_data(raw_data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('raw_data',@isstruct);
i_p.addParamValue('filter_set',NaN,@islogical);

i_p.parse(raw_data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use the filter_set variable to filter all the data sets before continuing
if (isempty(strmatch('filter_set',i_p.UsingDefaults)))
    these_names = fieldnames(raw_data);
    for j=1:size(these_names,1)
        raw_data.(these_names{j}) = raw_data.(these_names{j})(i_p.Results.filter_set,:);
    end
end

process_data = struct();

process_data.active_degrade = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.7;

process_data.active_degrade_10 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.1;
process_data.active_degrade_20 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.2;
process_data.active_degrade_30 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.3;
process_data.active_degrade_40 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.4;
process_data.active_degrade_50 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.5;
process_data.active_degrade_60 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.6;
process_data.active_degrade_70 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.7;
process_data.active_degrade_80 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.8;
process_data.active_degrade_90 = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
    & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < -0.9;

% process_data.active_degrade = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05 ...
%     & not(isnan(raw_data.cell_diff_percents)) & raw_data.cell_diff_percents < 0 ...
%     & not(isnan(raw_data.cell_diff_median)) & raw_data.cell_diff_median < 0 ...
%     & not(isnan(raw_data.cell_diff_mean)) & raw_data.cell_diff_mean < 0;

disp(['Detected ', num2str(sum(process_data.active_degrade(:))), ' invasion events.']);
% disp(['Bonferroni Corrected p-value threshold: ', num2str(0.05/bonferroni_correction)]);

process_data.live_cells = raw_data.tracking > -1;
process_data.longevities = sum(process_data.live_cells,2)/2;

process_data.ever_degrade = [];
for i=1:size(raw_data.tracking,1)
    process_data.ever_degrade = [process_data.ever_degrade, any(process_data.active_degrade(i,:))];
end

process_data.has_degraded = zeros(size(raw_data.tracking));
for i=1:size(raw_data.tracking,1)
    for j = 1:size(raw_data.tracking,2)
        process_data.has_degraded(i,j) = any(process_data.has_degraded(i,1:j));
    end
end

process_data.degrade_percentage = sum(process_data.has_degraded)/size(raw_data.tracking,1);
function find_full_exp_degrade_percents(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('min_longevity',10,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('..'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(exp_dir);
fields = filter_to_time_series(fields);

raw_data = struct('active_degrade',[],'longevity',[],'tracking',[]);

for j=1:length(fields)
    props_base = fullfile(exp_dir,fields(j).name,'cell_props');
    
    try
    degrade_data = csvread(fullfile(props_base,'active_degrade.csv'));
    raw_data.active_degrade = [raw_data.active_degrade; degrade_data];
    
    longev_data = csvread(fullfile(props_base,'longevity.csv'));
    raw_data.longevity = [raw_data.longevity;longev_data];
    
    tracking_data = csvread(fullfile(exp_dir,fields(j).name,'tracking_matrices','tracking_seq.csv'));
    raw_data.tracking = [raw_data.tracking;tracking_data];
    end
end

longev_filter = raw_data.longevity > i_p.Results.min_longevity;

processed_data = process_raw_data(raw_data,'filter_set',longev_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_dir = fullfile(exp_dir,'overall_results');
if (not(exist(output_dir,'dir')))
    mkdir(output_dir)
end
csvwrite(fullfile(output_dir,'degrade_percentage.csv'),processed_data.degrade_percentage);
csvwrite(fullfile(output_dir,'live_cells_per_image.csv'),processed_data.live_cells_per_image);

function processed_data = process_raw_data(raw_data,varargin)

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

if (isempty(strmatch('filter_set',i_p.UsingDefaults)))
    these_names = fieldnames(raw_data);
    for j=1:size(these_names,1)
        raw_data.(these_names{j}) = raw_data.(these_names{j})(i_p.Results.filter_set,:);
    end
end

processed_data = struct();

processed_data.live_cells = raw_data.tracking > -1;
processed_data.live_cells_per_image = sum(processed_data.live_cells);

processed_data.ever_degrade = [];
for i=1:size(raw_data.active_degrade,1)
    processed_data.ever_degrade = [processed_data.ever_degrade, any(raw_data.active_degrade(i,:))];
end

processed_data.has_degraded = zeros(size(raw_data.active_degrade));
for i=1:size(raw_data.active_degrade,1)
    for j = 1:size(raw_data.active_degrade,2)
        processed_data.has_degraded(i,j) = raw_data.active_degrade(i,j) | any(processed_data.has_degraded(i,1:j));
    end
end

processed_data.degrade_percentage = sum(processed_data.has_degraded)/size(raw_data.active_degrade,1);
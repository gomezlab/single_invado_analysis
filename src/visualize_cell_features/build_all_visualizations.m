function build_all_visualizations(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_reading = tic;
image_sets = cell(size(image_dirs,1),1);
for i = 1:size(image_dirs,1)
    image_sets{i} = read_in_file_set(fullfile(base_dir,image_dirs(i).name),filenames);
end
reading_time = toc(start_reading);
fprintf('It took %d seconds to read the image sets.\n',round(reading_time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Visualization Scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_tracking = tic;
try
    make_tracking_visualization(exp_dir,'image_sets',image_sets);
    fprintf('It took %d seconds to build the tracking vis.\n',round(toc(start_tracking)));
catch err
    fprintf('Problem with tracking: %s\n',err.message);
end

start_frames = tic;
try
    make_single_puncta_frames(exp_dir,'image_sets',image_sets);
    fprintf('It took %d seconds to build the single puncta frames.\n',round(toc(start_frames)));
catch err
    fprintf('Problem with single puncta frames: %s\n',err.message);
end

start_dual_highlight = tic;
try
    build_dual_highlight_image(exp_dir,'image_sets',image_sets);
    fprintf('It took %d seconds to build the side by side vis.\n',round(toc(start_dual_highlight)));
catch err
    fprintf('Problem with side by side visualizations: %s\n',err.message);
end

toc(start_time);
end

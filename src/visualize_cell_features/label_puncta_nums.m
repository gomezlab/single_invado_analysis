function label_puncta_nums(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('image_sets',NaN,@iscell);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'))
filenames = add_filenames_to_struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all the image directories
individual_images_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(individual_images_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

for i_num = 1:length(image_dirs)
    this_image_dir = fullfile(exp_dir,'individual_pictures',image_dirs(i_num).name);
    
    puncta_label = imread(fullfile(this_image_dir,'puncta_labeled.png'));
    
    size_increase_factor = 4;
    
    puncta_highlight = imread(fullfile(this_image_dir,'puncta_highlight.png'));
    puncta_highlight = imresize(puncta_highlight,size_increase_factor);
    imwrite(puncta_highlight,fullfile(this_image_dir,'puncta_nums.jpg'));

    gel_highlight = imread(fullfile(this_image_dir,'gel_highlight.png'));
    gel_highlight = imresize(gel_highlight,size_increase_factor);
    imwrite(gel_highlight,fullfile(this_image_dir,'gel_nums.jpg'));
    
    all_annotate = '';
    props = regionprops(puncta_label);
    for puncta_num = 1:max(puncta_label(:))
        pos_str = [' +',num2str(props(puncta_num).Centroid(1)*size_increase_factor),'+',num2str(props(puncta_num).Centroid(2)*size_increase_factor)];
        all_annotate = [all_annotate, ' -annotate', pos_str, ' ', num2str(puncta_num)]; %#ok<AGROW>
    end
    command_str = ['convert ', fullfile(this_image_dir,'puncta_nums.jpg'), ...
        ' -undercolor ''rgba(0,0,0,0.5)'' -font VeraBd.ttf -pointsize 16 -fill ''rgba(0,0,255,1)''', ...
        all_annotate, ' ', fullfile(this_image_dir,'puncta_nums.jpg'), ';'];

    system(command_str);
    
    command_str = ['convert ', fullfile(this_image_dir,'gel_nums.jpg'), ...
        ' -undercolor ''rgba(0,0,0,0.5)'' -font VeraBd.ttf -pointsize 16 -fill ''rgba(0,0,255,1)''', ...
        all_annotate, ' ', fullfile(this_image_dir,'gel_nums.jpg'), ';'];

    system(command_str);
end

toc(start_time);
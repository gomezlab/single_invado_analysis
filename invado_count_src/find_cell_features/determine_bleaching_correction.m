function determine_bleaching_correction(exp_dir,varargin)
% DETERMINE_BLEACHING_CORRECTION    searches through the gel images to
%                                   determine the photobleaching correction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;

i_p = inputParser;
i_p.FunctionName = 'DETERMINE_BLEACHING_CORRECTION';

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir);

i_p.addParamValue('output_dir',fullfile(exp_dir,'adhesion_props'),@ischar);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(exp_dir,varargin{:});

if (not(exist(i_p.Results.output_dir,'dir')))
    mkdir(i_p.Results.output_dir);
end

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine single image folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_dir = fullfile(exp_dir, 'individual_pictures');

single_image_folders = dir(image_dir);

assert(strcmp(single_image_folders(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(single_image_folders(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(single_image_folders(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

single_image_folders = single_image_folders(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the no cell region image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read in an image to get the size of the images to preallocate the size of
%the image that will keep track of which pixels have been in a cell or
%outside the registered range
test_image = imread(fullfile(image_dir,single_image_folders(1).name,filenames.cell_mask_filename));

no_cell_regions = ones(size(test_image));

for i=1:length(single_image_folders)
    cell_mask = imread(fullfile(image_dir,single_image_folders(i).name,filenames.cell_mask_filename));    
    no_cell_regions = no_cell_regions & not(cell_mask);
end

imwrite(no_cell_regions, fullfile(i_p.Results.output_dir,'no_cell_regions.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect the intensity correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gel_levels = zeros(length(single_image_folders),1);
gel_levels_outside_cell = zeros(size(gel_levels));

for i=1:length(single_image_folders)
    gel = imread(fullfile(image_dir,single_image_folders(i).name,filenames.gel_filename));
        
    gel_levels(i) = mean(gel(:));
    gel_levels_outside_cell(i) = mean(gel(no_cell_regions));    
end

for i=1:length(single_image_folders)
    dlmwrite(fullfile(image_dir, single_image_folders(i).name, filenames.intensity_correction_filename), ...
        1000/gel_levels_outside_cell(i));
end


%diagnostic plot
time_points = (0:(length(gel_levels) - 1))*5;
diag_fig_hnd = plot(time_points,gel_levels);
xlabel('Time (min)', 'Fontsize',16)
ylabel('Average Intensity', 'Fontsize',16);
hold on;
plot(time_points,gel_levels_outside_cell,'r');

y_limits = ylim();
ylim([0 y_limits(2)]);

legend('Overall','Outside Cell', 'location','SouthEast')
saveas(diag_fig_hnd,fullfile(i_p.Results.output_dir,'bleaching_curves.png'))
close all;

dlmwrite(fullfile(i_p.Results.output_dir,'bleaching_curves.csv'), ...
    [gel_levels_outside_cell,gel_levels]);

profile off;
if (i_p.Results.debug), profile viewer; end

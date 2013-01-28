function apply_bleaching_correction(exp_dir,varargin)
% APPLY_BLEACHING_CORRECTION    searches through the gel images and applies
%                               photobleaching correction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(i_p.Results.exp_dir,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%read in a cell mask to get the size of the images to preallocate the size of
%the matrix that will keep track of which pixels have been in a cell or
%outside the registered range
test_image = imread(fullfile(base_dir,image_dirs(1).name,filenames.gel));

no_cell_regions = ones(size(test_image));
for i=1:length(image_dirs)    
    if (exist(fullfile(base_dir,image_dirs(i).name,filenames.cell_mask),'file')) 
        cell_mask = imread(fullfile(base_dir,image_dirs(i).name,filenames.cell_mask));
        no_cell_regions = no_cell_regions & not(cell_mask);
    end
end
out_dir = fileparts(fullfile(base_dir,image_dirs(1).name,filenames.no_cell_regions));
if (not(exist(out_dir,'dir')))
    mkdir(out_dir);
end
imwrite(no_cell_regions, fullfile(base_dir,image_dirs(1).name,filenames.no_cell_regions));

%with the region of the image without cells found, we can then find the
%mean intensity outside the cells to determine an unbiased by degradation
%photobleaching correction
gel_levels_outside_cell = zeros(length(image_dirs),1);
gel_levels = zeros(length(image_dirs),1);
gel_levels_final = zeros(length(image_dirs),1);

for i=1:length(image_dirs)    
    gel = imread(fullfile(base_dir,image_dirs(i).name,filenames.gel));
    
    gel_levels(i) = mean(gel(:));
    gel_levels_outside_cell(i) = mean(gel(no_cell_regions));
    
    gel_corrected = gel .* (gel_levels_outside_cell(1)/gel_levels_outside_cell(i));
    gel_levels_final(i) = mean(gel_corrected(:));
    
    imwrite(gel_corrected,fullfile(base_dir,image_dirs(i).name,filenames.gel),'Bitdepth',16);
end

%diagnostic plot
time_points = (0:(length(gel_levels) - 1))*5;
diag_fig_hnd = plot(time_points,gel_levels);
xlabel('Time (min)', 'Fontsize',16)
ylabel('Average Intensity', 'Fontsize',16);
hold on;
plot(time_points,gel_levels_outside_cell,'r');
plot(time_points,gel_levels_final,'k');

y_limits = ylim();
ylim([0 y_limits(2)]);

legend('Overall','Outside Cell','Corrected','location','SouthEast')
print('-depsc2', fullfile(base_dir,image_dirs(1).name,filenames.bleaching_plot));
close all;

toc;
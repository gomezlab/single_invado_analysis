function determine_bleaching_correction(image_dir,varargin)
% DETERMINE_BLEACHING_CORRECTION    searches through the gel images to
%                                   determine the photobleaching correction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;

i_p = inputParser;
i_p.FunctionName = 'DETERMINE_BLEACHING_CORRECTION';

i_p.addRequired('image_dir',@(x)exist(x,'dir') == 7);

i_p.parse(image_dir);

i_p.addParamValue('output_dir',fullfile(image_dir,'..','adhesion_props'),@ischar);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(image_dir,varargin{:});

if (not(exist(i_p.Results.output_dir,'dir')))
    mkdir(i_p.Results.output_dir);
end

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull in data from the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_image_folders = dir(image_dir);

assert(strcmp(single_image_folders(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(single_image_folders(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(single_image_folders(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

single_image_folders = single_image_folders(3:end);

%read in an image to get the size of the images to preallocate the size of
%the image that will keep track of which pixels have been in a cell or
%outside the registered range
test_image = imread(fullfile(image_dir,single_image_folders(1).name,filenames.cell_mask_filename));

no_cell_regions = ones(size(test_image));
inside_registered = ones(size(test_image));

for i=1:length(single_image_folders)
    binary_shift = imread(fullfile(image_dir,single_image_folders(i).name,filenames.binary_shift_filename));
    cell_mask = imread(fullfile(image_dir,single_image_folders(i).name,filenames.cell_mask_filename));
    
    no_cell_regions = no_cell_regions & binary_shift;
    no_cell_regions = no_cell_regions & not(cell_mask);
    
    inside_registered = inside_registered & binary_shift;
end

imwrite(no_cell_regions, fullfile(i_p.Results.output_dir,'no_cell_regions.png'));

%check for the situation where there were less than 2% of pixels included
%in the no_cell_regions image, in that case, switch over to the image
%excluding the outside binary shift image
percent_outside = sum(sum(no_cell_regions))/(size(no_cell_regions,1)*size(no_cell_regions,2));

if (percent_outside < 0.02)
    no_cell_regions = inside_registered;
    imwrite(no_cell_regions, fullfile(i_p.Results.output_dir,'inside_registered.png'));
end

gel_levels = zeros(length(single_image_folders),1);
gel_levels_outside_cell = zeros(size(gel_levels));
gel_levels_puncta = zeros(size(gel_levels));

for i=1:length(single_image_folders)
    gel = imread(fullfile(image_dir,single_image_folders(i).name,filenames.gel_filename));
    scale_factor = double(intmax(class(gel)));
    gel  = double(gel)/scale_factor;
    
    puncta = imread(fullfile(image_dir,single_image_folders(i).name,filenames.puncta_filename));
    puncta = im2bw(puncta,0);
    
    gel_levels(i) = mean(gel(:));
    gel_levels_puncta(i) = mean(gel(puncta));
    gel_levels_outside_cell(i) = mean(gel(no_cell_regions));
    
    dlmwrite(fullfile(image_dir, single_image_folders(i).name, 'intensity_correction.csv'), ...
        0.1/mean(gel(no_cell_regions)));
end

%diagnostic plot
time_points = (0:(length(gel_levels) - 1))*5;
diag_fig_hnd = plot(time_points,gel_levels)
xlabel('Time (min)', 'Fontsize',16)
ylabel('Average Intensity', 'Fontsize',16);
hold on;
plot(time_points,gel_levels_puncta,'g');
plot(time_points,gel_levels_outside_cell,'r');

y_limits = ylim();
ylim([0 y_limits(2)]);

legend('Overall','Puncta','Outside Cell', 'location','SouthEast')
saveas(diag_fig_hnd,fullfile(i_p.Results.output_dir,'bleaching_curves.png'))
close all;

dlmwrite(fullfile(i_p.Results.output_dir,'bleaching_curves.csv'), ...
    [gel_levels_outside_cell,gel_levels, gel_levels_puncta]);

profile off;
if (i_p.Results.debug), profile viewer; end
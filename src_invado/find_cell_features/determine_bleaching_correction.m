function determine_bleaching_correction(image_dir,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_PUNCTA_PROPERTIES';

i_p.addRequired('image_dir',@(x)exist(x,'dir') == 7);

i_p.parse(image_dir);

i_p.addParamValue('output_dir',fullfile(image_dir,'..','adhesion_props'),@ischar);
i_p.addParamValue('adhesions_filename','puncta_labeled.png',@ischar);
i_p.addParamValue('puncta_filename','registered_focal_image.png',@ischar);
i_p.addParamValue('gel_filename','registered_gel.png',@ischar);
i_p.addParamValue('binary_shift_filename','binary_shift.png',@ischar);
i_p.addParamValue('cell_mask_filename','cell_mask.png',@ischar);

i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(image_dir,varargin{:});

if (not(exist(i_p.Results.output_dir,'dir')))
    mkdir(i_p.Results.output_dir);
end

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
test_image = imread(fullfile(image_dir,single_image_folders(1).name,i_p.Results.cell_mask_filename));

no_cell_regions = ones(size(test_image));
inside_registered = ones(size(test_image));

for i=1:length(single_image_folders)
    binary_shift = imread(fullfile(image_dir,single_image_folders(i).name,i_p.Results.binary_shift_filename));
    cell_mask = imread(fullfile(image_dir,single_image_folders(i).name,i_p.Results.cell_mask_filename));
    
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

gel_levels_outside_cell = zeros(length(single_image_folders),1);
gel_levels_puncta = zeros(length(single_image_folders),1);
gel_levels = zeros(length(single_image_folders),1);

for i=1:length(single_image_folders)
    gel = imread(fullfile(image_dir,single_image_folders(i).name,i_p.Results.gel_filename));
    scale_factor = double(intmax(class(gel)));
    gel  = double(gel)/scale_factor;
    
    puncta = imread(fullfile(image_dir,single_image_folders(i).name,i_p.Results.puncta_filename));
    puncta = im2bw(puncta,0);
    
    gel_levels(i) = mean(gel(:));
    gel_levels_puncta(i) = mean(gel(puncta));
    gel_levels_outside_cell(i) = mean(gel(no_cell_regions));
    
    dlmwrite(fullfile(image_dir, single_image_folders(i).name, 'intensity_correction.csv'), ...
        0.1/mean(gel(no_cell_regions)));
end

%diagnostic plot
% plot(gel_levels)
% xlabel('Time (min)')
% ylabel('Average Intensity');
% hold on;
% plot(gel_levels_puncta,'g');
% plot(gel_levels_outside_cell,'r');
% legend('Overall','Puncta','Outside Cell')

dlmwrite(fullfile(i_p.Results.output_dir,'bleaching_curves.csv'), ...
    [gel_levels; gel_levels_outside_cell; gel_levels_puncta]);
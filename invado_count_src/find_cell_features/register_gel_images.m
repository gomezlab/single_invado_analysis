function register_gel_images(I_file,reg_target_file,varargin)
% REGISTER_GEL_IMAGES    Performs a simple x-y translation registration on
%                        the two images provided, finding the optimum x-y
%                        shift by searching for the minimum in the mean
%                        difference
%
%   register_gel_images(I1,I2,OPTIONS) Register the first image (I1) to the
%   second image (I2)
%
%   Options:
%
%       -search_grid_resolution: number of pixels to jump in search for
%        optimum location on the -100 to 100 search grid space, defaults to
%        1
%       -output_dir: folder used to hold all the results, defaults to the
%        same folder as the image file, 'I1'
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile off; profile on;


i_p = inputParser;
i_p.FunctionName = 'REGISTER_GEL_IMAGES';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('reg_target_File',@(x)exist(x,'file') == 2);
i_p.addParamValue('search_grid_resolution',1, @(x)isnumeric(x) & x >= 1);
i_p.addParamValue('do_registration',1, @(x)isnumeric(x) && x == 1 || x == 0);
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,reg_target_file,varargin{:});

if (not(i_p.Results.do_registration))
    csvwrite(fullfile(i_p.Results.output_dir, 'affine_matrix.csv'), ...
        [cos(0) sin(0);-sin(0) cos(0); 0 0])
    return
end

%read in and normalize the input focal adhesion image
gel_image  = imread(I_file);
scale_factor = double(intmax(class(gel_image)));
gel_image  = double(gel_image)/scale_factor;

unreg_cell_mask = imread(fullfile(fileparts(I_file),'unreg_cell_mask.png'));

%read in and normalize the first image
reg_target  = imread(reg_target_file);
scale_factor = double(intmax(class(reg_target)));
reg_target  = double(reg_target)/scale_factor;

unreg_target_cell_mask = imread(fullfile(fileparts(reg_target_file),'unreg_cell_mask.png'));

%block out the regions around each of the cell masks
both_masks = unreg_cell_mask | unreg_target_cell_mask;
both_masks = imdilate(both_masks,strel('disk',40,0));

gel_image(both_masks) = NaN;
reg_target(both_masks) = NaN;

reg_target = reg_target .* (0.1/mean(reg_target(not(isnan(reg_target)))));
gel_image = gel_image .* (0.1/mean(gel_image(not(isnan(gel_image)))));

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row_shifts,col_shifts] = meshgrid(-50:i_p.Results.search_grid_resolution:50);
ms_diff = ones(size(row_shifts))*NaN;

%first we build out 5 layers of registration tests, 98% of the image
%registrations fall within this boundary
reg_layer_count = 0;
for temp = 1:5
    ms_diff = add_registration_layer(gel_image, reg_target, row_shifts, col_shifts, ms_diff);
    reg_layer_count = reg_layer_count + 1;
    if (i_p.Results.debug), disp(reg_layer_count); end
end

%if the best registration result is on the edge of the ms_diff matrix, we have
%to continue the search for the local minima, we don't want this search to
%continue forever though, so we cutoff the search when the ms_diff matrix is
%full, an unlikely but necessary check
while (best_reg_on_edge(ms_diff) && any(any(isnan(ms_diff))))
    ms_diff = add_registration_layer(gel_image, reg_target, row_shifts, col_shifts, ms_diff);
    reg_layer_count = reg_layer_count + 1;
    disp(reg_layer_count);
end

best_index = find(ms_diff == min(min(ms_diff)),1);
transform_matrix = [cos(0) sin(0);-sin(0) cos(0); row_shifts(best_index) col_shifts(best_index)];

% trimmed_ms_diff_binary = not(isnan(ms_diff));
% row_start = find(sum(trimmed_ms_diff_binary,1) > 0,1,'first');
% row_end = find(sum(trimmed_ms_diff_binary,1) > 0,1,'last');
% col_start = find(sum(trimmed_ms_diff_binary,2) > 0,1,'first');
% col_end = find(sum(trimmed_ms_diff_binary,2) > 0,1,'last');

trimmed_ms_diff = trim_nans(ms_diff);
assert(sum(sum(trimmed_ms_diff_binary)) == size(trimmed_ms_diff,1)*size(trimmed_ms_diff,2))
assert(all(all(isnumeric(trimmed_ms_diff))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csvwrite(fullfile(i_p.Results.output_dir, 'affine_matrix.csv'), transform_matrix)
csvwrite(fullfile(i_p.Results.output_dir, 'registration_diffs.csv'), trimmed_ms_diff)

if (i_p.Results.debug)
    profile viewer;
end
profile off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ms_diff = add_registration_layer(gel_image,reg_target,row_shifts,col_shifts, ms_diff)
% takes in the images, row/column shifts and the current ms_diff matrix and
% adds another layer of diff values

if (all(all(isnan(ms_diff))))
    middle_points_row = [ceil(size(row_shifts,1)/2), ceil(size(row_shifts,2)/2)];
    middle_points_col = [ceil(size(col_shifts,1)/2), ceil(size(col_shifts,2)/2)];
    
    this_row_shift = row_shifts(middle_points_row(1),middle_points_row(2));
    this_col_shift = col_shifts(middle_points_col(1),middle_points_col(2));
    
    transform = maketform('affine',[cos(0) sin(0);-sin(0) cos(0); this_row_shift this_col_shift]);
    
    binary_image = ones(size(gel_image));
    binary_shift = imtransform(binary_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
    gel_shift = imtransform(gel_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
    
    differences = reg_target.*binary_shift - gel_shift.*binary_shift;
    differences(~binary_shift) = NaN;
    differences = differences(not(isnan(differences)));
    ms_diff(middle_points_row(1),middle_points_row(2)) = mean(differences.*differences);
else
    %binary image showing which positions have been found
    ms_diff_calced = not(isnan(ms_diff));
    
    %expand the size of the binary image image by one pixel so that the
    %detected boundary only contains uncalculated diff points
    ms_diff_calced = bwperim(imdilate(ms_diff_calced, strel('square',3)));
    
    %scan through the shift matrices, finding the points to calc the new
    %diffs
    for i = 1:size(row_shifts,1)
        for j = 1:size(col_shifts,2)
            if (ms_diff_calced(i,j) == 0)
                continue;
            end
            
            this_row_shift = row_shifts(i,j);
            this_col_shift = col_shifts(i,j);
            
            transform = maketform('affine',[cos(0) sin(0);-sin(0) cos(0); this_col_shift this_row_shift]);
            
            binary_image = ones(size(gel_image));
            binary_shift = imtransform(binary_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
            gel_shift = imtransform(gel_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
            
            differences = reg_target - gel_shift;
            differences(~binary_shift) = NaN;
            differences = differences(not(isnan(differences)));
            ms_diff(i,j) = mean(differences.*differences);
        end
    end
end

function on_edge = best_reg_on_edge(ms_diff)
ms_diff_calced_border = bwperim(not(isnan(ms_diff)));

best_index = find(ms_diff == min(min(ms_diff)),1);

on_edge = ms_diff_calced_border(best_index);

function trimmed_ms_diff = trim_nans(ms_diff)

trimmed_ms_diff_binary = not(isnan(ms_diff));
row_start = find(sum(trimmed_ms_diff_binary,1) > 0,1,'first');
row_end = find(sum(trimmed_ms_diff_binary,1) > 0,1,'last');
col_start = find(sum(trimmed_ms_diff_binary,2) > 0,1,'first');
col_end = find(sum(trimmed_ms_diff_binary,2) > 0,1,'last');

trimmed_ms_diff = ms_diff(row_start:row_end, col_start:col_end);

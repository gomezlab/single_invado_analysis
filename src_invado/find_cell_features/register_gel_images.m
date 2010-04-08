function register_gel_images(I_file,reg_target,varargin)
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
profile on;
maxNumCompThreads(2);

i_p = inputParser;
i_p.FunctionName = 'REGISTER_GEL_IMAGES';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addRequired('reg_target',@(x)exist(x,'file') == 2);
i_p.addParamValue('search_grid_resolution',1, @(x)isnumeric(x) & x >= 1);
i_p.addParamValue('do_registration',1, @(x)isnumeric(x) & x == 1 || x == 0);
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,reg_target,varargin{:});

if (not(i_p.Results.do_registration))
    csvwrite(fullfile(i_p.Results.output_dir, 'affine_matrix.csv'), ... 
        [cos(0) sin(0);-sin(0) cos(0); 0 0])
    return
end

%read in and normalize the input focal adhesion image
gel_image  = imread(I_file);
scale_factor = double(intmax(class(gel_image)));
gel_image  = double(gel_image)/scale_factor;

%read in and normalize the first image
reg_target  = imread(reg_target);
scale_factor = double(intmax(class(reg_target)));
reg_target  = double(reg_target)/scale_factor;

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row_shifts,col_shifts] = meshgrid(-25:i_p.Results.search_grid_resolution:25);
ms_diff = zeros(size(row_shifts));

for i = 1:size(row_shifts,1)
    if (i_p.Results.debug)
        disp(col_shifts(i,1))
    end
    for j = 1:size(col_shifts,2)
        this_row_shift = row_shifts(i,j);
        this_col_shift = col_shifts(i,j);

        transform = maketform('affine',[cos(0) sin(0);-sin(0) cos(0); this_row_shift this_col_shift]);
        
        binary_image = ones(size(gel_image));
        binary_shift = imtransform(binary_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
        gel_shift = imtransform(gel_image, transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
        
        ms_diff(i,j) = sum(sum((reg_target.*binary_shift - gel_shift.*binary_shift).^2))/sum(sum(binary_shift));
    end
end
csvwrite(fullfile(i_p.Results.output_dir, 'registration_diffs.csv'), ms_diff)

best_index = find(ms_diff == min(min(ms_diff)),1);
transform_matrix = [cos(0) sin(0);-sin(0) cos(0); row_shifts(best_index) col_shifts(best_index)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csvwrite(fullfile(i_p.Results.output_dir, 'affine_matrix.csv'), transform_matrix)

if (i_p.Results.debug)
    profile viewer;
end
profile off;

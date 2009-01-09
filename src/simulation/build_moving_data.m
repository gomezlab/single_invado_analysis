function build_moving_data(varargin)
% BUILD_SIMULATED_DATA    Builds simulated focal adhesion data with known
%                         properties for use in analyzing the quality of
%                         the designed algorithms for studying focal
%                         adhesions
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'BUILD_MOVING_DATA';

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('output_dir', '.', @ischar);

i_p.parse(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ad_mean_intensity = 0.4237;

background_mean_intensity = 0.1220;

background_noise_var = 0.004;

max_ad_size = 10;

ad_padding = ceil(max_ad_size*0.8);

start_ad_size = 1;

ad_speed = 0.25;
image_width = ad_speed*25;

field_size = [max_ad_size^2 - max_ad_size*(start_ad_size - 1) + 2*ad_padding + ad_padding*(max_ad_size - start_ad_size + 1), 2*ad_padding + image_width + max_ad_size];

output_dir = fullfile('..','..','data','simulation',['moving_', num2str(ad_speed)],'Images','Paxillin');

%Create an image to be placed on the side apart from all the other images
side_ad_size = 3;
assert(side_ad_size <= max_ad_size, 'Side adhesion size must be less than or equal to the max adhesion size');
temp_ad = fspecial('gaussian',side_ad_size,side_ad_size/3);
temp_ad = temp_ad * (1/mean(mean(temp_ad))) * ad_mean_intensity;

side_image = zeros(field_size(1),2*ad_padding+size(temp_ad,2));
row_range = 1:size(temp_ad,1);
row_shift = floor((size(side_image,1) - size(temp_ad,1))/2);
row_range = row_range + row_shift;

col_range = 1:size(temp_ad,2);
col_shift = floor((size(side_image,2) - size(temp_ad,2))/2);
col_range = col_range + col_shift;

side_image(row_range,col_range) = temp_ad;

image_number = ceil(image_width/ad_speed)+2;

%Clear out the data directory
png_files = dir(fullfile(output_dir, '*.png'));
for i = 1:size(png_files,1), delete(fullfile(output_dir, png_files(i).name)); end

%Build the initial blank image
sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
imwrite([zeros(field_size), side_image], fullfile(output_dir,[sprintf(sprintf_format,1), '.png']))

for j = 2:(image_number-1)
    image = imnoise(zeros(field_size),'gaussian',background_mean_intensity,background_noise_var);
    for i=start_ad_size:max_ad_size
        temp_ad = fspecial('gaussian',i,i/3);
        temp_ad = temp_ad * (1/mean(mean(temp_ad))) * ad_mean_intensity;
        assert(size(temp_ad,1) == size(temp_ad,2));

        row_range = 1:size(temp_ad,1);
        row_shift = floor(ad_padding + ad_padding * (i - start_ad_size) + max_ad_size * (i - start_ad_size) + max_ad_size/2);
        row_range = row_range + row_shift;

        assert(size(row_range,2) == size(temp_ad,1), '%d, %d',size(row_range,2), size(temp_ad,1));
        assert(size(image,1) > row_range(end), 'Field size less than row range')
        
        col_range = 1:size(temp_ad,2);
        col_shift = floor(ad_padding + image_width * (j/image_number) + (max_ad_size - size(temp_ad,2))/2);
        col_range = col_range + col_shift;
        assert(size(col_range,2) == size(temp_ad,1), '%d, %d',size(col_range,2), size(temp_ad,1));
        assert(size(image,2) > col_range(end), 'Field size less than col range')
        
        image(row_range,col_range) = temp_ad;
    end

    sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
    if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
    imwrite([image, side_image], fullfile(output_dir,[sprintf(sprintf_format,j), '.png']))
end

sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
imwrite([zeros(field_size), side_image], fullfile(output_dir,[sprintf(sprintf_format,image_number), '.png']))
function build_stationary_data(varargin)
% BUILD_SIMULATED_DATA    Builds simulated focal adhesion data with known
%                         properties for use in analyzing the quality of
%                         the designed algorithms for studying focal
%                         adhesions
%
%   find_focal_adhesions(I,OPTIONS) locate the focal adhesions in image
%   file, 'I'
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'BUILD_SIMULATED_DATA';

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('output_dir', '.', @ischar);

i_p.parse(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('in_ad_pixels.mat')
% load('not_in_ad_pixels.mat')
% 
% in_ad_indexes = randperm(size(in_ad,1));
% not_in_ad_indexes = randperm(size(not_in_ad,1));


base_dir = '../../results/focal_adhesions/time_series_01/individual_pictures/';

% exp_dirs = dir('../../results/focal_adhesions/');
% assert(strmatch(exp_dirs(1).name,'.') && strmatch(exp_dirs(2).name,'..'))
% exp_dirs = exp_dirs(3:end);
% 
% for j = 1:length(exp_dirs)


dirs = dir(base_dir);
count = 0;

in_ad = [];
not_in_ad = [];
means = [];
for i=1:length(dirs)
    if (not(exist(fullfile(base_dir,dirs(i).name,'focal_image.png'),'file') && ...
            exist(fullfile(base_dir,dirs(i).name,'cell_mask.png'),'file') && ...
            exist(fullfile(base_dir,dirs(i).name,'adhesions_binary.png'),'file')))
        continue;
    end

    fi = imread(fullfile(base_dir,dirs(i).name,'focal_image.png'));
    fi = double(fi)/2^16;

    cm = imread(fullfile(base_dir,dirs(i).name,'cell_mask.png'));

    ads = imread(fullfile(base_dir,dirs(i).name,'adhesions_binary.png'));

    in_ad = [in_ad; fi(ads)];
    means(i-2,1) = mean(fi(ads));
    not_in_ad = [not_in_ad; fi(not(ads) & cm)];
    means(i-2,2) = mean(fi(not(ads) & cm));
end



ad_padding = 10;

max_ad_size = 10;

image_number = 100;

field_size = max_ad_size^2 + ad_padding*(max_ad_size+1);

binary_image = zeros(2*ad_padding+max_ad_size,field_size);

for i=1:max_ad_size
    temp_ad = zeros(max_ad_size);
    
    for j=1:size(temp_ad,1)
        for k=1:size(temp_ad,2)
            if (sqrt((j-size(temp_ad,1)/2)^2+(k-size(temp_ad,2)/2)^2) <= i/2)
                temp_ad(j,k) = 1;
            end
        end
    end

    assert(size(temp_ad,1) == size(temp_ad,2));

    row_range = floor(size(binary_image,1)/2)-floor((size(temp_ad,1) - 1)/2):floor(size(binary_image,1)/2)+ceil((size(temp_ad,1) - 1)/2);
    assert(size(row_range,2) == size(temp_ad,1), '%d, %d',size(row_range,2), size(temp_ad,1));

    col_range_start = ad_padding*i + max_ad_size*(i - 1) + (max_ad_size/2) - floor((size(temp_ad,1) - 1)/2);
    col_range = col_range_start : col_range_start + ceil((size(temp_ad,1) - 1)/2) + floor((size(temp_ad,1) - 1)/2);
    assert(size(col_range,2) == size(temp_ad,1), '%d, %d',size(col_range,2), size(temp_ad,1));

    binary_image(row_range,col_range) = temp_ad;
end

output_dir = fullfile('simulation','stationary','Images','Paxillin');

for i = 1:image_number
    output_image = zeros(size(binary_image));
    
    in_ad_pixels_to_fill = find(binary_image);
    
    for j = 1:size(in_ad_pixels_to_fill,1)
        try
            output_image(in_ad_pixels_to_fill(j)) = in_ad(in_ad_indexes(1));
            in_ad_indexes = in_ad_indexes(2:end);
        catch
            in_ad_indexes = randperm(size(in_ad,1));
            output_image(in_ad_pixels_to_fill(j)) = in_ad(in_ad_indexes(1));
            in_ad_indexes = in_ad_indexes(2:end);
        end
    end

    not_in_ad_pixels_to_fill = find(not(binary_image));

    for j = 1:size(not_in_ad_pixels_to_fill,1)
        try
            output_image(not_in_ad_pixels_to_fill(j)) = not_in_ad(not_in_ad_indexes(1));
            not_in_ad_indexes = not_in_ad_indexes(2:end);
        catch
            not_in_ad_indexes = randperm(size(not_in_ad,1));
            output_image(not_in_ad_pixels_to_fill(j)) = not_in_ad(not_in_ad_indexes(1));
            not_in_ad_indexes = not_in_ad_indexes(2:end);
        end
    end
    
%     for j=1:size(binary_image,1)
%         for k = 1:size(binary_image,2)
%             if (binary_image(j,k))
%                 try
%                     output_image(j,k) = in_ad(in_ad_indexes(1));
%                     in_ad = in_ad(2:end);
%                 catch
%                     in_ad_indexes = randperm(size(in_ad,1));
%                     output_image(j,k) = in_ad(in_ad_indexes(1));
%                     in_ad = in_ad(2:end);
%                 end
%             else
%                 try
%                     output_image(j,k) = not_in_ad(not_in_ad_indexes(1));
%                     not_in_ad = not_in_ad(2:end);
%                 catch
%                     in_ad_indexes = randperm(size(not_in_ad,1));
%                     output_image(j,k) = not_in_ad(not_in_ad_indexes(1));
%                     not_in_ad = not_in_ad(2:end);
%                 end
%             end
%         end
%         disp(j)
%     end
%     
    sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
    if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
    imwrite(output_image, fullfile(output_dir,[sprintf(sprintf_format,i), '.png']))
end


% ad_mean_intensity = 0.8;
% 
% background_mean_intensity = 0.1488;
% 
% noise_var = 0.004;
% 
% ad_padding = 10;
% 
% max_ad_size = 10;
% 
% field_size = max_ad_size^2 + ad_padding*(max_ad_size+1);
% 
% output_dir = fullfile('simulation','stationary','Images','Paxillin');
% 
% image_number = 100;
% 
% for j = 1:image_number
%     image = imnoise(zeros(2*ad_padding+max_ad_size,field_size),'gaussian',background_mean_intensity,noise_var);
%     for i=1:max_ad_size
%         temp_ad = fspecial('gaussian',i,i/3);
%         temp_ad = temp_ad * (1/max(max(temp_ad))) * 0.4384;
%         assert(size(temp_ad,1) == size(temp_ad,2));
% 
%         row_range = floor(size(image,1)/2)-floor((size(temp_ad,1) - 1)/2):floor(size(image,1)/2)+ceil((size(temp_ad,1) - 1)/2);
%         assert(size(row_range,2) == size(temp_ad,1), '%d, %d',size(row_range,2), size(temp_ad,1));
% 
%         col_range_start = ad_padding*i + max_ad_size*(i - 1) + (max_ad_size/2) - floor((size(temp_ad,1) - 1)/2);
%         col_range = col_range_start : col_range_start + ceil((size(temp_ad,1) - 1)/2) + floor((size(temp_ad,1) - 1)/2);
%         assert(size(col_range,2) == size(temp_ad,1), '%d, %d',size(col_range,2), size(temp_ad,1));
% 
%         image(row_range,col_range) = temp_ad;
%     end
%     
%     sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
%     if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
%     imwrite(image, fullfile(output_dir,[sprintf(sprintf_format,j), '.png']))
% end
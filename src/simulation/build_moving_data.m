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
i_p.addParamValue('output_dir', fullfile('..','..','data','simulation','moving_0_5','Images','Paxillin'), @ischar);

i_p.parse(varargin{:});

output_dir = i_p.Results.output_dir;

%Other Parameters

%source the 10th and 90th percentiles of all the adhesion intensities
min_ad_intensity = 0.2341950;
max_ad_intensity = 0.4723663;

speed = 0.5;

background_mean_intensity = 0.01;

background_noise_var = 0.005;

max_ad_size = 10;
min_ad_size = 1;

ad_padding = ceil(max_ad_size*0.4);

total_images = 25+2;
distance_covered = ceil((total_images - 2)*speed);

standard_frame_size = [max_ad_size + 2*ad_padding, max_ad_size + 2*ad_padding + distance_covered];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear out the directory currently holding the adhesions
png_files = dir(fullfile(output_dir, '*.png'));
for i = 1:size(png_files,1), delete(fullfile(output_dir, png_files(i).name)); end

%Create the side image that will appear in all the frames
side_image = zeros((max_ad_size-min_ad_size+1)*standard_frame_size(1), max_ad_size+4*ad_padding);
side_mid = [size(side_image)/2, size(side_image)/2];
if (mod(side_mid(1),1) ~= 0)
    side_mid(1) = side_mid(1) - 0.5;
end
if (mod(side_mid(2),1) ~= 0)
    side_mid(2) = side_mid(2) - 0.5;
end

this_ad = make_ad_matrix([floor(max_ad_size/2),floor(max_ad_size/2)], max_ad_intensity);
row_range = ceil(side_mid(1) - size(this_ad,1)/2):floor(side_mid(1) + size(this_ad,1)/2);
col_range = ceil(side_mid(2) - size(this_ad,2)/2):floor(side_mid(2) + size(this_ad,2)/2);
side_image(row_range, col_range) = this_ad;


for image_number = 2:(total_images - 1)
    %initialize all the individual adhesion frames
    image_frames = cell(max_ad_size-min_ad_size+1, 1);
    for i = 1:size(image_frames,1)
        for j = 1:size(image_frames,2)
            image_frames{i,j} = imnoise(zeros(standard_frame_size),'gaussian',background_mean_intensity,background_noise_var);
        end
    end
    
    %find the central point for this simulated adhesion
    ad_pos = [size(image_frames{i,j},1)/2, ... 
        (ad_padding + distance_covered*((image_number - 1)/(total_images - 2)) + max_ad_size/2)];
    if (mod(ad_pos(1),1) ~= 0)
        ad_pos(1) = ad_pos(1) - 0.5;
    end
    if (mod(ad_pos(2),1) ~= 0)
        ad_pos(2) = round(ad_pos(2));
    end
    
    %fill all the frames with adhesions
    for i = 1:size(image_frames,1)
        size_sequence = min_ad_size:max_ad_size;
        this_size = size_sequence(i);
        
        this_ad = make_ad_matrix([this_size,this_size],mean([max_ad_intensity,min_ad_intensity]));
        
        row_range = ceil(ad_pos(1) - size(this_ad,1)/2):floor(ad_pos(1) + size(this_ad,1)/2);
        if (mod(size(this_ad,1),2) == 0)
            row_range = row_range(1:(length(row_range) - 1));
        end
        col_range = ceil(ad_pos(2) - size(this_ad,2)/2):floor(ad_pos(2) + size(this_ad,2)/2);
        if (mod(size(this_ad,2),2) == 0)
            col_range = col_range(1:(length(col_range) - 1));
        end
        
        assert(size(this_ad,1) == length(row_range))
        assert(size(this_ad,2) == length(col_range));
        
        image_frames{i,j}(row_range,col_range) = image_frames{i,j}(row_range,col_range) + this_ad;
    end
    
    final_image = put_together_cell_images(image_frames);
    sprintf_format = ['%0', num2str(length(num2str(total_images))), 'd'];
    if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
    imwrite([final_image, side_image], fullfile(output_dir,[sprintf(sprintf_format,image_number), '.png']))
end

%Build the initial and final blank image
sprintf_format = ['%0', num2str(length(num2str(total_images))), 'd'];
if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
imwrite([zeros(size(put_together_cell_images(image_frames))), side_image], ... 
    fullfile(output_dir,[sprintf(sprintf_format,1), '.png']))
imwrite([zeros(size(put_together_cell_images(image_frames))), side_image], ... 
    fullfile(output_dir,[sprintf(sprintf_format,total_images), '.png']))



function overall_image = put_together_cell_images(image_set)

overall_image = [];
for i = 1:size(image_set,1)
    this_col = image_set{i,1};
    for j = 2:size(image_set,2)
        this_col = [this_col, image_set{i,j}]; %#ok<AGROW>
    end
    overall_image = [overall_image; this_col]; %#ok<AGROW>
end
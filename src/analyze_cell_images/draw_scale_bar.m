function image_with_bar = draw_scale_bar(varargin)

%%Parameter Checking and Setup
if (size(varargin,2) < 1)
    error('ERROR: draw_scale_bar - not enough parameters specified, at least one needed');
elseif (isstruct(varargin{1}))
    image_data = varargin{1};
    if (isfield(image_data,'pixel_size'))
        pixel_size = image_data.pixel_size;
    else
        error('ERROR: draw_scale_bar - the pixel size must be available');
    end

    if (isfield(image_data,varargin{2}))
        image_no_bar = getfield(image_data,varargin{2});
    else
        error('ERROR: draw_scale_bar - Can''t find image to put bar in');
    end

    if(isfield(image_data,'cell_mask'))
        cell_mask = image_data.cell_mask;
    else
        error('ERROR: draw_scale_bar - Cell mask needed');
    end

    bar_size = 10E-6;
    if (size(varargin,2) >= 3)
        if (isnumeric(varargin{3}))
            bar_size = varargin{3};
        else
            error('ERROR: draw_scale_bar - Expected third parameter to be a number, bar_size');
        end
    end

    position_code = 0;
    if (size(varargin,2) >= 4)
        if (isnumeric(varargin{4}))
            position_code = varargin{4};
        else
            error('ERROR: draw_scale_bar - Expected fourth parameter to be a number, position_code');
        end
    end
end

%%Global Variables
offset_standard = 0.03;

%%Main Program
[image_height image_width image_dim] = size(image_no_bar);

bar_width = round(bar_size/pixel_size);
bar_height = round(0.01*image_height);
bar = ones(bar_height,bar_width);
bar_size = size(bar);

%%%%%Determine Position Code
best_position_code = 2;

if (not(position_code))
    position_code = best_position_code;
end

switch position_code
    case 1
        row_offset = offset_standard;
        col_offset = 1-offset_standard;
    case 2
        row_offset = 1-offset_standard;
        col_offset = 1-offset_standard;
    case 3
        row_offset = 1-offset_standard;
        col_offset = offset_standard;
    case 4
        row_offset = offset_standard;
        col_offset = offset_standard;
    otherwise
        error('ERROR: draw_scale_bar - unexpected position code, must be between 0-4');
end

%%%%%% Draw the image with the bar
image_with_bar = image_no_bar;

for i=1:image_dim
    switch position_code
        case 1
            bar_rows = round(row_offset*image_height):round(row_offset*image_height)-1+bar_size(1);
            bar_cols = round(col_offset*image_width)-bar_size(2):round(col_offset*image_width)-1;
        case 2
            bar_rows = round(row_offset*image_height)-bar_size(1):round(row_offset*image_height)-1;
            bar_cols = round(col_offset*image_width)-bar_size(2):round(col_offset*image_width)-1;
        case 3
            bar_rows = round(row_offset*image_height)-bar_size(1):round(row_offset*image_height)-1;
            bar_cols = round(col_offset*image_width):round(col_offset*image_width)-1+bar_size(2);
        case 4
            bar_rows = round(row_offset*image_height):round(row_offset*image_height)-1+bar_size(1);
            bar_cols = round(col_offset*image_width):round(col_offset*image_width)-1+bar_size(2);
        otherwise
            error('ERROR: draw_scale_bar - unexpected position code, must be between 0-4');
    end
    image_with_bar(bar_rows,bar_cols,i) = bar;
end

end
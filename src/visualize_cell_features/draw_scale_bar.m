function image_with_bar = draw_scale_bar(image_no_bar,bar_size)
%DRAW_SCALE_BAR    Adds a scale bar to an image
%
%   Position Codes:
%       -1: Upper right hand corner
%       -2: Lower right hand corner
%       -3: Lower left hand corner
%       -4: Upper left hand corner
%   

%%Parameter Checking and Setup

position_code = 2;
pixel_size = 0.215051;

%%Global Variables
offset_standard = 0.03;

%%Main Program
[image_height image_width image_dim] = size(image_no_bar);

bar_width = round(bar_size/pixel_size);
bar_height = round(0.01*image_height);
bar = 255*ones(bar_height,bar_width);
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
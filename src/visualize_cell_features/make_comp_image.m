function composite_image = make_comp_image(varargin);

%%focal_edge_highlights,focal_image,cell_mask
%%Error Checking and Setup
if (isstruct(varargin{1}))
    image_data = varargin{1};
    if (isfield(image_data,'focal_edge_highlights'))
        focal_edge_highlights = image_data.focal_edge_highlights;
    else
        error('ERROR: make_comp_image - focal_edge_highlights struct entry not found');
    end

    if (isfield(image_data,'focal_image'))
        focal_image = image_data.focal_image;
    else
        error('ERROR: make_comp_image - focal_image struct entry not found');
    end

    if (isfield(image_data,'cell_mask'))
        cell_mask = image_data.cell_mask;
    else
        error('ERROR: make_comp_image - cell_mask struct entry not found');
    end
elseif (isnumeric(varargin{1}))
    focal_edge_highlights = varargin{1};
    if (size(varargin,2) >= 3)
        if (isnumeric(varargin{2}))
            focal_image = varargin{2};
        else
            error('ERROR: make_comp_image - must provide focal_image as second parameter when first is matrix');
        end
        
        if (isnumeric(varargin{3}) || islogical(varargin{3}))
            cell_mask = varargin{3};
        else
            error('ERROR: make_comp_image - must provide cell_mask as third parameter when first is matrix');
        end
            
    else
        error('ERROR: make_comp_image - must provide three parameters when first is a matrix');
    end
else
    error('ERROR: make_comp_image - provide appropriate struct or numeric matrices');
end

%%Main Program
col_sums = sum(cell_mask);
row_sums = sum(cell_mask');

min_row = find(row_sums,1,'first') - 10;
max_row = find(row_sums,1,'last') + 10;
min_col = find(col_sums,1,'first') - 10;
max_col = find(col_sums,1,'last') + 10;

if ((min_row - 10) < 0)
    min_row = 1;
end

if ((max_row + 10) > size(cell_mask,1))
    max_row = size(cell_mask,1);
end

if ((min_col - 10) < 0)
    min_col = 1;
end

if ((max_col + 10) > size(cell_mask,2))
    max_col = size(cell_mask,2);
end

if (length(size(focal_image)) < 3)
    focal_image = cat(3,focal_image,focal_image,focal_image);
end

if (length(size(focal_edge_highlights)) < 3)
    focal_edge_highlights = cat(3,focal_edge_highlights,focal_edge_highlights,focal_edge_highlights);
end

clipped_highlights = focal_edge_highlights(min_row:max_row,min_col:max_col,1:3);
clipped_focal = focal_image(min_row:max_row,min_col:max_col,1:3);
spacer = ones(size(clipped_focal,1),floor(0.03*size(clipped_focal,2)),3);

composite_image = [clipped_highlights, spacer, clipped_focal];

end
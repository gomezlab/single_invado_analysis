function composite_image = make_comp_image(focal_edge_highlights,focal_image,cell_mask);

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
    temp = focal_image;
    focal_image = ones(size(focal_image,1),size(focal_image,2),3);
    focal_image(:,:,1) = temp;
    focal_image(:,:,2) = temp;
    focal_image(:,:,3) = temp;
end

clipped_highlights = focal_edge_highlights(min_row:max_row,min_col:max_col,1:3);
clipped_focal = focal_image(min_row:max_row,min_col:max_col,1:3);
spacer = ones(size(clipped_focal,1),floor(0.03*size(clipped_focal,2)),3);

composite_image = [clipped_highlights, spacer, clipped_focal];

end
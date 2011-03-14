function trimmed_image = remove_region_outside_registered(image,binary_shift)

min_row = find(sum(binary_shift,2),1,'first');
max_row = find(sum(binary_shift,2),1,'last');
min_col = find(sum(binary_shift),1,'first');
max_col = find(sum(binary_shift),1,'last');

trimmed_image = image(min_row:max_row, min_col:max_col);
function [trimmed_image] = trim_to_edge(I,Mask)

col_sums = sum(Mask);
row_sums = sum(Mask');

min_row = find(row_sums,1,'first') - 10;
max_row = find(row_sums,1,'last') + 10;
min_col = find(col_sums,1,'first') - 10;
max_col = find(col_sums,1,'last') + 10;

if ((min_row - 10) < 0)
    min_row = 1;
end

if ((max_row + 10) > size(Mask,1))
    max_row = size(Mask,1);
end

if ((min_col - 10) < 0)
    min_col = 1;
end

if ((max_col + 10) > size(Mask,2))
    max_col = size(Mask,2);
end

if (length(size(I)) < 3)
    I = cat(3,I,I,I);
end

trimmed_image = I(min_row:max_row,min_col:max_col,1:3);
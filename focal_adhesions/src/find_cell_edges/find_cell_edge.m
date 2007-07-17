function [edge_binary_image, counts, x] = find_cell_edge(I,range)

I = imfilter(I,fspecial('gaussian',[15 15],1));

edge_binary_image = zeros(size(I,1),size(I,2));

sorted_pixel_vals = sort(I(:));

sorted_pixel_vals = sorted_pixel_vals(floor(0.005*length(I(:))):length(sorted_pixel_vals));

[counts,x] = hist(sorted_pixel_vals,256);

smoothed_counts = smooth(counts);

max_index = find(counts == max(counts))
if (size(max_index) > 1) 
    max_index = max_index(1);
end

line(x(1:length(x)),counts(1:length(x)));
hold on; line(x(1:length(x)),smoothed_counts(1:length(x)),'color','r');

count_diffs = zeros(1,length(max_index:length(counts)));
length(count_diffs)
for i = max_index:(length(counts) - 1)
    count_diffs(i) = counts(i) - counts(i+1);
end
count_diffs = count_diffs/max(count_diffs(:));



end
base_folder = '../../data/time_points/';

edge_image = imread(strcat(base_folder,'1/N-myr mRFP.tif'),3);
%edge_image = edge_image(1:1040,401:1392);
%edge_image = edge_image(1:1040,1:599);
%edge_image = edge_image(201:1040,1:599);
%edge_image = edge_image(1:600,1:599);
edge_range = double(max(edge_image(:)) - min(edge_image(:)));

focal_image = imread(strcat(base_folder,'1/EGFP-Paxillin.tif'),3);

edge_image = normalize_grayscale_images(edge_image);
focal_image = normalize_grayscale_images(focal_image);

%[binary_edge,counts,x] = find_cell_edge(edge_image,edge_range);

I = edge_image;
range = edge_range;

I = imfilter(I,fspecial('gaussian',[15 15],1));

sorted_pixel_vals = sort(I(:));

sorted_pixel_vals = sorted_pixel_vals(floor(0.005*length(I(:))):length(sorted_pixel_vals));

[counts,x] = hist(sorted_pixel_vals,10000);

smoothed_counts = smooth(counts);

max_index = find(counts == max(counts));
if (size(max_index) > 1) 
    max_index = max_index(1);
end

min_count_val = min(smoothed_counts(max_index:floor(max_index+0.2*range)));
min_count_val_index = find(smoothed_counts < max(smoothed_counts)*0.001,1,'first');
if (size(min_count_val_index,1) > 1) 
    min_count_val_index = min_count_val_index(1);
end

threshold_val = x(min_count_val_index);
threshold_val = 0.08;

edge_binary_image = bwperim(im2bw(I,threshold_val));

binary_edge = edge_binary_image;

temp = zeros(size(edge_image,1),size(edge_image,2),3);
temp(:,:,1) = edge_image;
temp(:,:,2) = edge_image;
temp(:,:,3) = edge_image;

highlight_focal = zeros(size(focal_image,1),size(focal_image,2),3);
highlight_focal(:,:,1) = focal_image;
highlight_focal(:,:,2) = focal_image;
highlight_focal(:,:,3) = focal_image;

for i=1:size(edge_image,1)
    for j=1:size(edge_image,2)
        if (binary_edge(i,j))
            temp(i,j,2) = 1;
            highlight_focal(i,j,2) = 1;
        end
    end
end
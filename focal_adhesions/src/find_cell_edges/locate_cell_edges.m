base_folder = '../../data/time_points/';

edge_image = imread(strcat(base_folder,'1/N-myr mRFP.tif'),1);
edge_range = double(max(edge_image(:)) - min(edge_image(:)));

focal_image = imread(strcat(base_folder,'1/EGFP-Paxillin.tif'),1);

edge_image = normalize_grayscale_images(edge_image);
focal_image = normalize_grayscale_images(focal_image);

threshold_val = adaptive_thresh(I)

edge_binary_image = bwperim(im2bw(I,threshold_val));

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
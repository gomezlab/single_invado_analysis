base_folder = '../../data/time_points/';

edge_image = imread(strcat(base_folder,'1/N-myr mRFP.tif'),1);
focal_image = imread(strcat(base_folder,'1/EGFP-Paxillin.tif'),1);

edge_image = normalize_grayscale_images(edge_image);
focal_image = normalize_grayscale_images(focal_image);

threshold_val = adaptive_thresh(edge_image);

edge_binary_image = bwperim(im2bw(edge_image,threshold_val));

highlighted_edge_image = create_highlighted_image(edge_image,edge_binary_image);
highlighted_focal_image = create_highlighted_image(focal_image,edge_binary_image);

imwrite(highlighted_edge_image,strcat(base_folder,'edge.png'));
imwrite(highlighted_focal_image,strcat(base_folder,'focal.png'));
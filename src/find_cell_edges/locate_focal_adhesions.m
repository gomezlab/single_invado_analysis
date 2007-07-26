base_folder = '../../data/time_points/';

number_of_timepoints = 5;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end
    
    image_set_cell_number = size(imfinfo(strcat(base_folder,num2str(i),'/N-myr mRFP.tif')),2);
    %image_set_cell_number = 04;
    
    for j = 1:image_set_cell_number
    %for j = 4:image_set_cell_number
        padded_cell_num = sprintf(strcat('%0', num2str(length(num2str(image_set_cell_number))), 'd'),j);
        padded_time_point_num = sprintf(strcat('%0', num2str(length(num2str(number_of_timepoints))), 'd'),i);
        output_directory = strcat(base_folder,'each_cell/',padded_time_point_num,'/',padded_cell_num,'/');

        focal_image = imread(strcat(base_folder,num2str(i),'/EGFP-Paxillin.tif'),j);
        focal_image = normalize_grayscale_image(focal_image);
        
        if (exist(strcat(output_directory,'cell_mask.png')))
            cell_mask = imread(strcat(output_directory,'cell_mask.png'));
        else
            edge_image = normalize_grayscale_image(imread(strcat(base_folder,num2str(i),'/N-myr mRFP.tif'),j));
            edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));            
            edge_binary_image = clean_up_edge_image(edge_binary_image);
            cell_mask = imfill(edge_binary_image,'holes');
        end
        
        focal_image = focal_image.*cell_mask;
        
        se = strel('disk',2);
        Ie = imerode(focal_image, se);
        Iobr = imreconstruct(Ie, focal_image);
        
        Iobrd = imdilate(Iobr, se);
        Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
        Iobrcbr = imcomplement(Iobrcbr);
        
        fgm = imregionalmax(Iobrcbr);
        
        inverted_focal_image = -focal_image + 1;
        inverted_focal_image = imfilter(inverted_focal_image,fspecial('gaussian'));
        inverted_focal_image(~cell_mask) = 1.1;
        inverted_focal_image = imimposemin(inverted_focal_image,fgm);
        
        watershed_labels = watershed(inverted_focal_image);
        watershed_labels(~cell_mask) = 0;
        watershed_labels(bwperim(cell_mask)) = 0;
        
        watershed_edges = zeros(size(focal_image,1),size(focal_image,2));
        watershed_edges(find(watershed_labels >= 1)) = 0;
        watershed_edges(find(watershed_labels == 0)) = 1;
        watershed_edges(~cell_mask) = 0;
        focal_edge_highlights = create_highlighted_image(focal_image,watershed_edges);
        
        full_cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
        full_cell_centroid = [full_cell_centroid.Centroid];
        full_cell_centroid = round(full_cell_centroid);
        focal_edge_highlights(full_cell_centroid(2)-5:full_cell_centroid(2)+5,full_cell_centroid(1)-5:full_cell_centroid(1)+5,1) = 1;
        
        centroid_stats = regionprops(watershed_labels,'centroid');
        centroid_stats = [centroid_stats.Centroid];
        
        dists = bwdist(~cell_mask);
        max_dist = max(dists(:));
        min_dist = min(dists(:));
        
        for k = 1:2:size(centroid_stats,2)
            middle_point_col = round(centroid_stats(k));
            middle_point_row = round(centroid_stats(k+1));
            col_range = middle_point_col-1:middle_point_col+1;
            row_range = middle_point_row-1:middle_point_row+1;
            
            this_dist = dists(middle_point_row,middle_point_col)/max_dist;
            
            color = ones(3,3,3);
            color(:,:,1) = 1-this_dist;
            color(:,:,2) = zeros(3,3);
            color(:,:,3) = this_dist;
            
            focal_edge_highlights(row_range,col_range,1:3) = color;
        end
        
        if (not(exist(output_directory,'dir')))
            mkdir(output_directory);
        end
        
        imwrite(focal_edge_highlights,strcat(output_directory,'focal_edges.png'));

        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end
    end
end
profile off;
run_information = profile('info');

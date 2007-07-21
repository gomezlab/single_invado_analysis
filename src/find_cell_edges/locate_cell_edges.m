base_folder = '../../data/time_points/';

number_of_cells = 19;
number_of_timepoints = 5;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end
    for j = 1:number_of_cells
        image_num = j;
        edge_image = imread(strcat(base_folder,num2str(i),'/N-myr mRFP.tif'),j);
        focal_image = imread(strcat(base_folder,num2str(i),'/EGFP-Paxillin.tif'),j);

        edge_image = normalize_grayscale_image(edge_image);
        focal_image = normalize_grayscale_image(focal_image);

        edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));
        
        original_binary_image = edge_binary_image;
        
        edge_binary_image = clean_up_edge_image(edge_binary_image);

        highlighted_edge_image = create_highlighted_image(edge_image,edge_binary_image);
        highlighted_original_image = create_highlighted_image(edge_image,original_binary_image);
        side_by_side_edge_image = [highlighted_edge_image, ones(size(edge_image,1),floor(0.05*size(edge_image,2)),3), highlighted_original_image];
        
        highlighted_focal_image = create_highlighted_image(focal_image,edge_binary_image);

        padded_cell_num = sprintf(strcat('%0', num2str(length(num2str(number_of_cells))), 'd'),j);
        output_directory = strcat(base_folder,'each_cell/',padded_cell_num,'/');
        if (not(exist(output_directory,'dir')))
            mkdir(output_directory);
        end
        
        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d',j)
            end
        end
        
        padded_time_point_num = sprintf(strcat('%0', num2str(length(num2str(number_of_timepoints))), 'd'),i);
        imwrite(highlighted_edge_image,strcat(output_directory,'edge',padded_time_point_num,'.png'));
        imwrite(side_by_side_edge_image,strcat(output_directory,'side_by_side',padded_time_point_num,'.png'));
        imwrite(highlighted_focal_image,strcat(output_directory,'focal',padded_time_point_num,'.png'));
        imwrite(imfill(edge_binary_image,'holes'),strcat(output_directory,'cell_mask',padded_time_point_num,'.png'));
    end
end
profile off;
p = profile('info');
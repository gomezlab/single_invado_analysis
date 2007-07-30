base_folder = '../../data/time_points/';

number_of_timepoints = 5;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end
    
    image_set_cell_number = size(imfinfo([base_folder,num2str(i),'/N-myr mRFP.tif']),2);
    
    for j = 1:image_set_cell_number
        edge_image = imread([base_folder,num2str(i),'/N-myr mRFP.tif'],j);
        focal_image = imread([base_folder,num2str(i),'/EGFP-Paxillin.tif'],j);

        edge_image = normalize_grayscale_image(edge_image);
        focal_image = normalize_grayscale_image(focal_image);

        edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));
        
        original_binary_image = edge_binary_image;
        
        edge_binary_image = clean_up_edge_image(edge_binary_image);

        highlighted_edge_image = create_highlighted_image(edge_image,edge_binary_image);
        highlighted_original_image = create_highlighted_image(edge_image,original_binary_image);
        side_by_side_edge_image = [highlighted_edge_image, ones(size(edge_image,1),floor(0.05*size(edge_image,2)),3), highlighted_original_image];     
        highlighted_focal_image = create_highlighted_image(focal_image,edge_binary_image);

        padded_cell_num = sprintf(['%0', num2str(length(num2str(image_set_cell_number))), 'd'],j);
        padded_time_point_num = sprintf(['%0', num2str(length(num2str(number_of_timepoints))), 'd'],i);
        output_directory = [base_folder,'individual_pictures/',padded_time_point_num,'/',padded_cell_num,'/'];
        if (not(exist(output_directory,'dir')))
            mkdir(output_directory);
        end

        imwrite(highlighted_edge_image,[output_directory,'edge.png']);
        imwrite(side_by_side_edge_image,[output_directory,'side_by_side.png']);
        imwrite(highlighted_focal_image,[output_directory,'focal.png']);
        imwrite(imfill(edge_binary_image,'holes'),[output_directory,'cell_mask.png']);
        
        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end
    end
end
profile off;
run_information = profile('info');
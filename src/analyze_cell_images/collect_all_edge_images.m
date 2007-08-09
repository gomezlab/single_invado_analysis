base_folder = '../../data/time_points/';

number_of_timepoints = 1;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end

    image_set_cell_number = size(imfinfo([base_folder,num2str(i),'/N-myr mRFP.tif']),2);

    for j = 1:image_set_cell_number
        image_data.edge_image = imread([base_folder,num2str(i),'/N-myr mRFP.tif'],j);
        image_data.focal_image = imread([base_folder,num2str(i),'/EGFP-Paxillin.tif'],j);

        image_data.edge_image = normalize_grayscale_image(image_data.edge_image);
        image_data.focal_image = normalize_grayscale_image(image_data.focal_image);

        image_data.edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(image_data.edge_image,0.2)));

        image_data.edge_binary_image = clean_up_edge_image(image_data.edge_binary_image);

        image_data.highlighted_edge_image = create_highlighted_image(image_data.edge_image,image_data.edge_binary_image);
        image_data.highlighted_focal_image = create_highlighted_image(image_data.focal_image,image_data.edge_binary_image);

        image_data.padded_cell_num = sprintf(['%0', num2str(length(num2str(image_set_cell_number))), 'd'],j);
        image_data.padded_time_point_num = sprintf(['%0', num2str(length(num2str(number_of_timepoints))), 'd'],i);
        image_data.output_directory = [base_folder,'individual_pictures/',image_data.padded_time_point_num,'/',image_data.padded_cell_num,'/'];
        if (not(exist(image_data.output_directory,'dir')))
            mkdir(image_data.output_directory);
        end

        imwrite(image_data.highlighted_edge_image,[image_data.output_directory,'edge.png']);
        imwrite(image_data.highlighted_focal_image,[image_data.output_directory,'focal.png']);
        imwrite(imfill(image_data.edge_binary_image,'holes'),[image_data.output_directory,'cell_mask.png']);

        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end
    end
end
profile off;
run_information = profile('info');
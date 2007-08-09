base_folder = fullfile('..','..','data','time_points');

number_of_timepoints = 5;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end

    image_set_cell_number = size(imfinfo(fullfile(base_folder,num2str(i),'N-myr mRFP.tif')),2);

    for j = 1:image_set_cell_number
        image_data.padded_cell_num = sprintf(['%0', num2str(length(num2str(image_set_cell_number))), 'd'],j);
        image_data.padded_time_point_num = sprintf(['%0', num2str(length(num2str(number_of_timepoints))), 'd'],i);
        image_data.output_directory = fullfile(base_folder,'individual_pictures',image_data.padded_time_point_num,image_data.padded_cell_num);
        
        image_data.edge_image = imread(fullfile(base_folder,num2str(i),'/N-myr mRFP.tif'),j);
        image_data.focal_image = imread(fullfile(base_folder,num2str(i),'/EGFP-Paxillin.tif'),j);

        image_data.edge_image = normalize_grayscale_image(image_data.edge_image);
        image_data.focal_image = normalize_grayscale_image(image_data.focal_image);

        image_data.cell_mask = create_cell_edge_image(image_data.edge_image,image_data.output_directory);
        image_data.edge_binary_image = bwperim(image_data.cell_mask);
        
        image_data.highlighted_edge_image = create_highlighted_image(image_data.edge_image,image_data.edge_binary_image);
        image_data.highlighted_focal_image = create_highlighted_image(image_data.focal_image,image_data.edge_binary_image);

        if (not(exist(image_data.output_directory,'dir')))
            mkdir(image_data.output_directory);
        end

        imwrite(image_data.highlighted_edge_image,fullfile(image_data.output_directory,'edge.png'));
        imwrite(image_data.highlighted_focal_image,fullfile(image_data.output_directory,'focal.png'));

        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end
    end
end
profile off;
run_information = profile('info');
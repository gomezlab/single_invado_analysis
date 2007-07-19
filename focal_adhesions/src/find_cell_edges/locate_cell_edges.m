base_folder = '../../data/time_points/';

number_of_cells = 19;
number_of_timepoints = 5;
debug = 1;

for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end
    for j = 1:number_of_cells
        image_num = j;
        edge_image = imread(strcat(base_folder,num2str(i),'/N-myr mRFP.tif'),image_num);
        focal_image = imread(strcat(base_folder,num2str(i),'/EGFP-Paxillin.tif'),image_num);

        edge_image = normalize_grayscale_images(edge_image);
        focal_image = normalize_grayscale_images(focal_image);

        edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image)));
        
        edge_binary_image = clean_up_edge_image(edge_binary_image);

        highlighted_edge_image = create_highlighted_image(edge_image,edge_binary_image);
        highlighted_focal_image = create_highlighted_image(focal_image,edge_binary_image);

        padded_cell_num = sprintf(strcat('%0', num2str(length(num2str(number_of_cells))), 'd'),j);
        output_directory = strcat(base_folder,'each_cell/',padded_cell_num,'/');
        if (not(exist(output_directory)))
            mkdir(output_directory);
        end
        
        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d',j)
            end
        end
        
        padded_time_point_num = sprintf(strcat('%0', num2str(length(num2str(number_of_timepoints))), 'd'),i);
        imwrite(highlighted_edge_image,strcat(output_directory,'edge',padded_time_point_num,'.png'));
        imwrite(highlighted_focal_image,strcat(output_directory,'focal',padded_time_point_num,'.png'));
    end
end
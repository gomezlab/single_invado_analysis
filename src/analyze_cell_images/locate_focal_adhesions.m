base_folder = '../../data/time_points/';

number_of_timepoints = 5;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end

    image_set_cell_number = size(imfinfo(strcat(base_folder,num2str(i),'/N-myr mRFP.tif')),2);
    %image_set_cell_number = 4;

    for j = 1:image_set_cell_number
    %for j = 4:image_set_cell_number
        
        padded_cell_num = sprintf(['%0', num2str(length(num2str(image_set_cell_number))), 'd'],j);
        padded_time_point_num = sprintf(['%0', num2str(length(num2str(number_of_timepoints))), 'd'],i);
        output_directory = [base_folder,'individual_pictures/',padded_time_point_num,'/',padded_cell_num,'/'];

        focal_image = normalize_grayscale_image(imread([base_folder,num2str(i),'/EGFP-Paxillin.tif'],j));

        if (exist(strcat(output_directory,'cell_mask.png')))
            cell_mask = imread(strcat(output_directory,'cell_mask.png'));
        else
            edge_image = normalize_grayscale_image(imread([base_folder,num2str(i),'/N-myr mRFP.tif'],j));
            edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));
            edge_binary_image = clean_up_edge_image(edge_binary_image);
            cell_mask = imfill(edge_binary_image,'holes');
        end

        focal_image = focal_image.*cell_mask;

        focal_markers = find_focal_adhesion_markers(focal_image,cell_mask);

        inverted_focal_image = -focal_image + 1;
        inverted_focal_image = imfilter(inverted_focal_image,fspecial('gaussian'));
        inverted_focal_image(~cell_mask) = 1.1;
        inverted_focal_image = imimposemin(inverted_focal_image,focal_markers);

        watershed_labels = watershed(inverted_focal_image);
        watershed_labels(~cell_mask) = 0;
        watershed_labels(bwperim(cell_mask)) = 0;

        watershed_edges = zeros(size(focal_image,1),size(focal_image,2));
        watershed_edges(find(watershed_labels >= 1)) = 0;
        watershed_edges(find(watershed_labels == 0)) = 1;
        watershed_edges(~cell_mask) = 0;
        focal_edge_highlights = create_highlighted_image(focal_image,watershed_edges);

        focal_edge_highlights = draw_centroid_dots(focal_edge_highlights, cell_mask, watershed_labels);
        [focal_edge_highlights,identified_adhesions] = find_each_watershed_adhesion(focal_image, focal_edge_highlights, watershed_labels, cell_mask);

        labeled_adhesions = bwlabel(identified_adhesions);
        
        adhesion_props = regionprops(labeled_adhesions,'all');
        
        if (not(exist(output_directory,'dir')))
            mkdir(output_directory);
        end
        
        if (not(exist([base_folder,'all'],'dir')))
            mkdir([base_folder,'all']);
        end
        
        composite_image = make_comp_image(focal_edge_highlights,focal_image,cell_mask);
        
        %imwrite(focal_edge_highlights,[output_directory,'focal_edges.png']);
        imwrite(focal_edge_highlights,[base_folder,'all/','focal_edges_',padded_time_point_num,'_',padded_cell_num,'.png']);
        %imwrite(composite_image,[output_directory,'comp.png']);

        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end
    end
end
profile off;
run_information = profile('info');

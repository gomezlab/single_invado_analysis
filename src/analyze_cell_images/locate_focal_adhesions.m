base_folder = '../../data/time_points/';

number_of_timepoints = 1;
debug = 1;
profile on;
for i = 1:number_of_timepoints
    if (debug)
        sprintf('Time Point Number: %02d',i)
    end

    image_set_cell_number = size(imfinfo(strcat(base_folder,num2str(i),'/N-myr mRFP.tif')),2);
    image_set_cell_number = 4;

    %for j = 1:image_set_cell_number
    for j = 4:image_set_cell_number
        
        image_data.pixel_size = importdata([base_folder,'pixel_size.txt']);
        image_data.padded_cell_num = sprintf(['%0', num2str(length(num2str(image_set_cell_number))), 'd'],j);
        image_data.padded_time_point_num = sprintf(['%0', num2str(length(num2str(number_of_timepoints))), 'd'],i);
        image_data.output_directory = [base_folder,'individual_pictures/',image_data.padded_time_point_num,'/',image_data.padded_cell_num,'/'];
        image_data.interesting_directory = [base_folder,'many_adhesions/',image_data.padded_time_point_num,'/',image_data.padded_cell_num,'/'];

        image_data.original_focal_image = normalize_grayscale_image(imread([base_folder,num2str(i),'/EGFP-Paxillin.tif'],j));

        if (exist(strcat(image_data.output_directory,'cell_mask.png')))
            image_data.cell_mask = imread(strcat(image_data.output_directory,'cell_mask.png'));
            image_data.cell_edge = bwperim(image_data.cell_mask);
        else
            edge_image = normalize_grayscale_image(imread([base_folder,num2str(i),'/N-myr mRFP.tif'],j));
            edge_binary_image = bwperim(im2bw(edge_image,adaptive_thresh(edge_image,0.2)));
            edge_binary_image = clean_up_edge_image(edge_binary_image);
            image_data.cell_mask = imfill(edge_binary_image,'holes');
            image_data.cell_edge = bwperim(image_data.cell_mask);
        end

        image_data.focal_image = image_data.original_focal_image.*image_data.cell_mask;

        image_data.focal_markers = find_focal_adhesion_markers(image_data);

        image_data.inverted_focal_image = -image_data.focal_image + 1;
        image_data.inverted_focal_image = imfilter(image_data.inverted_focal_image,fspecial('gaussian'));
        image_data.inverted_focal_image(~image_data.cell_mask) = 1.1;
        image_data.inverted_focal_image = imimposemin(image_data.inverted_focal_image,image_data.focal_markers);

        image_data.watershed_labels = watershed(image_data.inverted_focal_image);
        image_data.watershed_labels(~image_data.cell_mask) = 0;
        image_data.watershed_labels(image_data.cell_edge) = 0;

        image_data.watershed_edges = zeros(size(image_data.focal_image,1),size(image_data.focal_image,2));
        image_data.watershed_edges(find(image_data.watershed_labels >= 1)) = 0;
        image_data.watershed_edges(find(image_data.watershed_labels == 0)) = 1;
        image_data.watershed_edges(~image_data.cell_mask) = 0;
        
        image_data.identified_adhesions = find_watershed_adhesions(image_data);
        
        %image_data.focal_edge_highlights = create_highlighted_image(image_data,'watershed_edges');
        %image_data.focal_edge_highlights = draw_centroid_dots(image_data);
        image_data.focal_edge_highlights = create_highlighted_image(image_data,'identified_adhesions',1);
        image_data.focal_edge_highlights = draw_scale_bar(image_data,'focal_edge_highlights');

        image_data.labeled_adhesions = bwlabel(image_data.identified_adhesions);
        image_data.adhesion_props = collect_adhesion_properties(image_data);

        image_data.composite_image = make_comp_image(image_data);
        image_data.composite_image = draw_scale_bar(image_data,'composite_image');

        if (not(exist(image_data.output_directory,'dir')))
            mkdir(image_data.output_directory);
        end
        imwrite(image_data.focal_edge_highlights,[image_data.output_directory,'focal_highlights.png']);
        imwrite(image_data.composite_image,[image_data.output_directory,'comp.png']);
        imwrite(image_data.focal_edge_highlights,[image_data.output_directory,'all_adhesions.png']);

        if (max(image_data.labeled_adhesions(:)) > 300)
            if (not(exist(image_data.interesting_directory,'dir')))
                mkdir(image_data.interesting_directory);
            end
            imwrite(image_data.focal_edge_highlights,[image_data.interesting_directory,'focal_highlights.png']);
            imwrite(image_data.composite_image,[image_data.interesting_directory,'comp.png']);
        end

        if (debug)
            if (mod(j,5) == 0)
                sprintf('Cell Number: %02d / %02d',j,image_set_cell_number)
            end
        end

        %clear image_data;
    end
end
profile off;
run_information = profile('info');

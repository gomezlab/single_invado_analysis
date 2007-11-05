function make_movie(cfg_file)

[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

%this_movie = avifile(fullfile(movie_path,movie_filename));

addpath(genpath(path_folders));

exclude_fullfile = fullfile(fileparts(original_i_file),excluded_frames_filename);
if (exist(exclude_fullfile,'file'))
    excluded_frames = load(exclude_fullfile);
else
    excluded_frames = 0;
end

if (exist('out_path','var'))
    if (not(exist(out_path,'dir')))
        mkdir(out_path);
    end
end

extr_val_fullfile = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_fullfile,'file'))
    min_max = load(extr_val_fullfile);
end

i_count = size(imfinfo(original_i_file),2);

c_map = jet(i_count);

i_size = size(imread(original_i_file,1));
pix_count = i_size(1)*i_size(2);
edge_image = ones(i_size(1),i_size(2),3);
edge_image_ad = ones(i_size(1),i_size(2),3);

test_seq = load('test_seq.csv') + 1;
i_seen = 0;
for i = 1:i_count
    if (find(i==excluded_frames))
        continue;
    end
    i_seen = i_seen + 1;
    
    orig_i = normalize_grayscale_image(imread(original_i_file,i),min_max(1),min_max(2));

    i_size = size(orig_i);

    padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);

    I_1 = imread(fullfile(I_folder_1,padded_i_num,filename_1));
    
    I_2 = imread(fullfile(I_folder_2,padded_i_num,filename_2));
    markers = imread(fullfile(I_folder_1,padded_i_num,'focal_markers.png'));

    I_lab = bwlabel(I_2);
    I_test = zeros(i_size(1),i_size(2));
    I_test(find(I_lab == test_seq(i_seen))) = 1;
    highlighted_test = create_highlighted_image(orig_i,bwperim(I_test));
    highlighted_test = create_highlighted_image(highlighted_test,bwperim(I_2 & not(I_test)),1);    

    highlighted_2 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(highlighted_1,bwperim(I_1 & not(I_2)),1);
    %highlighted_1 = create_highlighted_image(highlighted_1,bwperim(markers),3);

    cell_edge_1 = bwperim(imread(fullfile(I_folder_1,padded_i_num,edge_filename)));

    for j = 1:3
        edge_image(find(cell_edge_1)+(j-1)*pix_count) = c_map(i,j);
        edge_image_ad(find(cell_edge_1)+(j-1)*pix_count) = c_map(i,j);
        edge_image_ad(find(bwperim(I_2))+(j-1)*pix_count) = c_map(i,j);
    end

    %frame = [highlighted_1,0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_2];
    frame = [edge_image_ad,0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_test];

    if (exist('out_path','var'))
        imwrite(frame,fullfile(out_path,[out_prefix,padded_i_num,'.png']));
    end

    %this_movie = addframe(this_movie,frame);
end

system('/Users/mbergins/bin/ffmpeg -r 5 -b 8000 -y -i ../../results/norm/time_series_1/movie_files/time%03d.png ../../results/norm/time_series_1/test.mp4');
function make_movie_frames(cfg_file)

[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

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

all_seqs = zeros(281,198);
seq_map = jet(281);
for i = 0:280
    all_seqs(i + 1,:) = load(fullfile('track_seqs',['seq_',num2str(i),'.csv'])) + 1;
end
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
    I_test(find(I_lab == all_seqs(1,i_seen))) = 1;
    highlighted_test = create_highlighted_image(orig_i,bwperim(I_test));
    highlighted_test = create_highlighted_image(highlighted_test,bwperim(I_2 & not(I_test)),'color',1);    

    highlighted_all = create_highlighted_image(orig_i,bwperim(I_2));
    for j=1:281
        this_high = zeros(i_size(1),i_size(2));
        this_high(find(I_lab == all_seqs(j,i_seen))) = 1;
        highlighted_all = create_highlighted_image(highlighted_all,bwperim(this_high),'color',seq_map(j,:));
    end
    
    highlighted_2 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(orig_i,bwperim(I_2));
    highlighted_1 = create_highlighted_image(highlighted_1,bwperim(I_1 & not(I_2)),'color',1);
    %highlighted_1 = create_highlighted_image(highlighted_1,bwperim(markers),3);

    cell_edge_1 = bwperim(imread(fullfile(I_folder_1,padded_i_num,edge_filename)));

    for j = 1:3
        edge_image(find(cell_edge_1)+(j-1)*pix_count) = c_map(i,j);
        edge_image_ad(find(cell_edge_1)+(j-1)*pix_count) = c_map(i,j);
        edge_image_ad(find(bwperim(I_2))+(j-1)*pix_count) = c_map(i,j);
    end

    %frame = [highlighted_1,0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_2];
    frame = [cat(3,orig_i,orig_i,orig_i),0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_all];

    if (exist('out_path','var'))
        output_filename = fullfile(out_path,[out_prefix,padded_i_num,'.png']);
        fullpath = fileparts(output_filename);
        if (not(exist(fullpath,'dir')))
            mkdir(fullpath);
        end
        imwrite(frame,output_filename);
    end

end

%system('/Users/mbergins/bin/ffmpeg -sameq -r 5 -b 8000 -y -i ../../results/norm/time_series_1/movie_files/time%03d.png -f huffyuv ../../results/norm/time_series_1/test.mp4');
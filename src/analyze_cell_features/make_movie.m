function make_movie(cfg_file)

[cfg_file_path,cfg_filename] = fileparts(cfg_file);
addpath(cfg_file_path);
eval(cfg_filename);
rmpath(cfg_file_path);

this_movie = avifile(fullfile(movie_path,movie_filename));

addpath(genpath(path_folders));

exclude_fullfile = fullfile(fileparts(original_i_file),excluded_frames_filename);
if (exist(exclude_fullfile,'file'))
    excluded_frames = load(exclude_fullfile);
else
    excluded_frames = 0;
end

extr_val_fullfile = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_fullfile,'file'))
    min_max = load(extr_val_fullfile);
end

i_count = size(imfinfo(original_i_file),2);

c_map = jet(i_count);

for i = 1:i_count
   if (find(i==excluded_frames))
       continue;
   end
   orig_i = normalize_grayscale_image(imread(original_i_file,i),min_max(1),min_max(2));
   
   padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
   
   I_1 = imread(fullfile(I_folder_1,padded_i_num,filename_1));
   
   
   I_2 = imread(fullfile(I_folder_2,padded_i_num,filename_2));
   markers = imread(fullfile(I_folder_1,padded_i_num,'focal_markers.png'));

   highlighted_2 = create_highlighted_image(orig_i,bwperim(I_2));
   highlighted_1 = create_highlighted_image(orig_i,bwperim(I_2));
   highlighted_1 = create_highlighted_image(highlighted_1,bwperim(I_1 & not(I_2)),1);
   %highlighted_1 = create_highlighted_image(highlighted_1,bwperim(markers),3);
   
   frame = [highlighted_1,0.5*ones(size(orig_i,1),round(0.05*size(orig_i,2)),3),highlighted_2];
   this_movie = addframe(this_movie,frame);
end

this_movie = close(this_movie);
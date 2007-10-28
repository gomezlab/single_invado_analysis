original_i_file = '../../data/time_series_1/Paxillin.tif';

I_folder_1 = '../../results/norm/time_series_1/individual_pictures/';
filename_1 = 'adhesions.png';
I_folder_2 = '../../results/norm/time_series_1/individual_pictures/';
filename_2 = 'time_filtered.png';

excluded_frames_filename = 'excluded_frames.txt';
extr_val_filename        = 'min_max_vals.txt';

movie_filename = 'both.avi';
movie_path = '.';

path_folders = '..';


addpath(genpath(path_folders));

exclude_fullfile = fullfile(fileparts(original_i_file),excluded_frames_filename);
if (exist(exclude_fullfile,'file'))
    excluded_frames = load(exclude_fullfile)
end

extr_val_fullfile = fullfile(fileparts(original_i_file),extr_val_filename);
if (exist(extr_val_fullfile,'file'))
    min_max = load(extr_val_fullfile);
end

i_count = size(imfinfo(original_i_file),2);

c_map = jet(i_count);

for i = 1:i_count
   if (find(i==excluded_frames))
       i
       continue;
   end
   orig_i = normalize_grayscale_image(imread(original_i_file),min_max(1),min_max(2));
   
   padded_i_num = sprintf(['%0',num2str(length(num2str(i_count))),'d'],i);
   
   I_1 = imread(fullfile(I_folder_1,padded_i_num,filename_1));
   I_2 = imread(fullfile(I_folder_2,padded_i_num,filename_2));
end
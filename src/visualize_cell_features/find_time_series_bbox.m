function bbox = find_time_series_bbox(folder,varargin)
% FIND_BINARY_BOUNDING_BOX    finds the minimum box which completely encloses
%                             the provided binary image
%
%   [min_x,min_y,max_x,max_y] = find_binary_bounding_box(I) finds the min
%   and max, x and y coordinates that completely enclose the binary 
%   image, 'I'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('folder',@(x)exist(x,'dir')==7);
i_p.addParamValue('image_filename','cell_mask.png',@ischar);

i_p.parse(folder,varargin{:});

image_filename = i_p.Results.image_filename;

folders_to_exclude = {'.','..'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_files = dir(folder);
num_files = size(all_files,1);

bbox = [Inf, Inf, -Inf, -Inf];
min_x = Inf;
min_y = Inf;
max_x = -Inf;
max_y = -Inf;

for i = 1:num_files
    full_folder_path = fullfile(folder,all_files(i).name);
    if (not(exist(full_folder_path,'dir')))
        continue;
    end
    if (strmatch(all_files(i).name,folders_to_exclude))
        continue;
    end
    
    image_path = fullfile(full_folder_path,image_filename);
    
    if (not(exist(image_path,'file')))
        continue;
    end
    
    this_bbox = find_binary_bounding_box(imread(image_path));
    
    if (bbox(1) > this_bbox(1)), bbox(1) = this_bbox(1); end
    if (bbox(2) > this_bbox(2)), bbox(2) = this_bbox(2); end
    if (bbox(3) < this_bbox(3)), bbox(3) = this_bbox(3); end
    if (bbox(4) < this_bbox(4)), bbox(4) = this_bbox(4); end
end
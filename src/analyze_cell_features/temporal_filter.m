function temporal_filter(folder,num,degree,varargin)
% TEMPORAL_FILTER    filters a given binary feature image based on the
%                    presence of overlapping features in a set of images
%
%   Note: This script assumes 
%
%   temporal_filter(F_1,F_2,Options) given two files, F_1 and F_2, which contain
%   binary image features, find the features in F_1 that overlap the
%   features in F_2, output the filtered features in the same folder as
%   F_1, but with the file name 'time_filtered.png', the additional options
%   in 'Options' will also be used
%
%   Options:
%       -'out_dir': specifies a directory where the filtered image will be
%        output, short option 'd', has priority over 'd'
%
%       -'out_filename': specifies a filename to use for the filtered
%        output image, short option 'f', has priority over 'f'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'TEMPORAL_FILTER';

i_p.addRequired('folder',@(x)exist(x,'dir')==7);
i_p.addRequired('num',@isnumeric);
i_p.addRequired('degree',@isnumeric);

i_p.parse(folder,num,degree);

folder_guess = fullfile(folder,num2str(num));
leading_zeros = 0;
while exist(folder_guess,'dir') ~= 7
    leading_zeros = leading_zeros + 1;
    if (leading_zeros > 50)
        break;
    end
    folder_guess = fullfile(folder,sprintf(['%0',num2str(leading_zeros),'d'],num));
end

i_p.addOptional('out_dir',folder_guess,@(x)exist(x,'dir')==7);
i_p.addOptional('out_filename','time_filtered.png',@ischar);
i_p.addOptional('in_filename','adhesions.png',@ischar);
i_p.parse(folder,num,degree,varargin{:});

%collect all the parsed input values into shorter variable names and read
%in the feature images

top_folder = i_p.Results.folder;
num = i_p.Results.num;
degree = i_p.Results.degree;
out_dir = i_p.Results.out_dir;
out_file = i_p.Results.out_filename;
in_file = i_p.Results.in_filename;

base_image_file = fullfile(top_folder,num2str(num),in_file);
leading_zeros = 0;
while exist(base_image_file,'file') ~= 2
    leading_zeros = leading_zeros + 1;
    if (leading_zeros > 50)
        error('TEMPORAL_FILTER: can''t find base image file');
    end
    base_image_file = fullfile(top_folder,sprintf(['%0',num2str(leading_zeros),'d'],num),in_file);
end

padded_base_num = sprintf(['%0',num2str(leading_zeros),'d'],num);
base_image = imread(base_image_file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_set = ones(size(base_image,1),size(base_image,2),degree);

offset = 0;
for i=1:degree
    if (exist(fullfile(top_folder,sprintf('%03d',num+i+offset),'adhesions.png'),'file'))
        image_set(:,:,i) = imread(fullfile(top_folder,sprintf('%03d',num+i+offset),'adhesions.png'));
    end 
end

intersected_image = base_image & image_set(:,:,1);
for i=2:degree
    intersected_image = intersected_image & image_set(:,:,i);
end

labeled_base = bwlabel(base_image);
out_image = zeros(size(base_image));
c = 0;

for i=1:max(labeled_base(:))
    if (sum(intersected_image(find(labeled_base==i))))
        out_image(find(labeled_base==i)) = 1;
    else
        c = c + 1;
    end
end
disp(c)

imwrite(out_image,fullfile(out_dir,out_file));

end
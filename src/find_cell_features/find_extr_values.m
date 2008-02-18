function min_max = find_extr_values(varargin)
% FIND_EXTR_VALUES     find the maximum and minumum values in a provided
%                      image or set of stacked tiff images
%
%   [min_max] = find_extr_values(I_stack_file) finds the maximum and
%   minumum value in all of the images in the file 'I_file', which is
%   either a single image or a stack of tiff images, a two element array is
%   returned, the first element in the minimum and the second is the
%   maximum
%
%   [min_max] = find_extr_values(I_file_1,I_file_2,...) finds the maximum
%   and minumum value in all of the images, 'I_file_1, I_file_2, etc',
%   included on the command line, a two element array is returned, the
%   first element in the minimum and the second is the maximum
%
%   [min_max] = find_extr_values(Folder,I_file_1,I_file_2,...) finds the
%   maximum and minumum value in all of the images, 'I_file_1, I_file_2,
%   etc', included on the command line, a two element array is returned,
%   the first element is the minimum and the second is the maximum, the
%   minimum and maximum are also written to 'min_max_vals.csv' in the
%   folder, 'Folder', in csv format 
%
%   [min_max] = find_extr_values(Folder,I_stack_file) finds the maximum and
%   minumum value in all of the images in the file 'I_file', which is
%   either a single image or a stack of tiff images, a two element array is
%   returned, the first element in the minimum and the second is the
%   maximum, the minimum and maximum are also written to 'min_max_vals.csv'
%   in the folder, 'Folder', in csv format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_output_file = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_max = [Inf, 0];

if (exist(fileparts(varargin{1}),'dir'))
    output_file = varargin{1};
    write_output_file = 1;
end

for i = 1:size(varargin,2)
    if (write_output_file && i == 1)
        continue;
    end
    
    file = varargin{i};
    
    for j=1:size(imfinfo(varargin{i}),2)
        temp = imread(file,j);

        if (min(temp(:)) < min_max(1))
            min_max(1) = min(temp(:));
        end
        if (max(temp(:)) > min_max(2))
            min_max(2) = max(temp(:));
        end
    end
end

min_max = double(min_max);

if (write_output_file)
    csvwrite(output_file,min_max);
end
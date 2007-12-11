function min_max = find_extr_values(file)
% FIND_EXTR_VALUES     find the maximum and minumum values in a provided
%                      image or set of stacked tiff images 
%
%   [min_max] = find_extr_values(I_stack_file) finds the maximum and
%   minumum value in all of the images in the file 'I_file', which is
%   either a single image or a stack of tiff images, a two element arrary is 
%   returned, the first element in the minimum and the second is the maximum  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_EXTR_VALS';

i_p.addRequired('file',@(x)exist(x,'file')==2);

i_p.parse(file);

stack_file = i_p.Results.file;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_max = [Inf, 0];
for i=1:size(imfinfo(stack_file),2)
    if (size(imfinfo(stack_file),2) == 1)
        temp = imread(stack_file);
    else
        temp = imread(stack_file,i);
    end

    if (min(temp(:)) < min_max(1))
        min_max(1) = min(temp(:));
    end
    if (max(temp(:)) > min_max(2))
        min_max(2) = max(temp(:));
    end
end

min_max = double(min_max);

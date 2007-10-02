function [min_val max_val] = find_extr_values(file)
% FIND_MAX_VALUE_IN_STACK     find the maximum and minumum values in a
%                             provided image or set of stacked tiff images
%
%   [min,max] = find_max_value_in_stack(I_stack_file) finds the maximum and
%   minumum value in all of the images in the file 'I_file', which is
%   either a single image or a stack of tiff images


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'find_extr_values_in_stack';
i_p.addRequired('file',@ischar);

i_p.parse(file);

stack_file = i_p.Results.file;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_val = 0;
min_val = inf;
for i=1:size(imfinfo(stack_file),2)
    if (size(imfinfo(stack_file),2) == 1)
        temp = imread(stack_file);
    else
        temp = imread(stack_file,i);
    end
    
    if (max(temp(:)) > max_val)
        max_val = max(temp(:));
    end
    if (min(temp(:)) < min_val)
        min_val = min(temp(:));
    end
end

min_val = double(min_val);
max_val = double(max_val);
function write_normalized_image(I_file,out_file,varargin)
% WRITE_GRAYSCALE_IMAGE   Write a normalized image to a provided output
%                         file using the minimum and maximum values
%                         specified in a provided file
%
%   write_grayscale_image(I,out) normalize the image file 'I', writing the 
%   result to the file 'out'
%
%   write_grayscale_image(I,out'I_num',num) normalize image number 'num' in the 
%   stacked image file 'I', using the value in the file 'min_max', writing the 
%   result to the file 'out'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'WRITE_NORMALIZED_IMAGE';

i_p.addRequired('I_file',@(x)exist(x,'file')==2);

i_p.addRequired('out_file',@(x)exist(fileparts(x),'dir')==7);

i_p.addParamValue('I_num',1,@(x)isnumeric(x) && x>0);
i_p.addParamValue('debug',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));
i_p.addParamValue('d',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));

i_p.parse(I_file,out_file,varargin{:});

%Determine if the image file specified has more than one image embedded, if
%so, make sure there is a 'I_num' parameter
if (size(imfinfo(I_file),2) > 1)
    if (any(strcmp('I_num',i_p.UsingDefaults)))
        error(['ERROR: ',i_p.FunctionName,' - Image file specified has more than one image embedded, must specify an ''I_num'' parameter']);
    end
    I_num = i_p.Results.I_num;
else
    I_num = 1;
end

input_image = imread(I_file,I_num);

output_bits = 8;
if (intmax(class(input_image)) > 255) 
    output_bits = 16;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

possible_ranges = [2^8, 2^12, 2^16];
max_val = possible_ranges(find(possible_ranges > max(input_image(:)), 1, 'first')); 

normalized_image = normalize_grayscale_image(input_image,'min_max',[0,max_val - 1]);

imwrite(normalized_image,out_file,'Bitdepth',output_bits);

end

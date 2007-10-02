function [I] = normalize_grayscale_image(varargin)
% NORMALIZE_GRAYSCALE_IMAGE    using a given image, the same image with the
%                              minimum and maximum pixel value adjusted to
%                              0 and 1 respectively is returned in double
%                              format
%
%   normalize_grayscale_image(I) normalize the image 'I' to have min and
%   max values of 0 and 1
%
%   normalize_grayscale_image(I,max,min) normalize the image 'I' to have
%   min and max values of 'min' and 'max'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isempty(varargin))
    error('ERROR: normalize_grayscale_image - command line options must be provided, see ''help normalize_grayscale_image''');
elseif (size(varargin,2) == 1)
    if (not(isnumeric(varargin{1})))
        error('ERROR: normalize_grayscale_image - first command line options must be numeric, see ''help normalize_grayscale_image''');
    end
    I = double(varargin{1});
    min_val = min(I(:));
    max_val = max(I(:))-min_val;
elseif (size(varargin,2) == 3)
    if (not(isnumeric(varargin{1})))
        error('ERROR: normalize_grayscale_image - first command line options must be numeric, see ''help normalize_grayscale_image''');
    end
    if (not(isnumeric(varargin{2})))
        error('ERROR: normalize_grayscale_image - second command line options must be numeric, see ''help normalize_grayscale_image''');
    end
    if (not(isnumeric(varargin{3})))
        error('ERROR: normalize_grayscale_image - third command line options must be numeric, see ''help normalize_grayscale_image''');
    end    
    
    I = double(varargin{1});
    min_val = varargin{2};
    max_val = varargin{3};
else
    error('ERROR: normalize_grayscale_image - unexpected number of command line options provided, see ''help normalize_grayscale_image''');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size = size(I);

if (size(image_size,2) > 2)
    if (image_size(3) > 1)
        I = I(:,:,1);
    end
end

I = I - min_val;
I = I/max_val;

end
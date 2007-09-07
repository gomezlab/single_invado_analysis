function threshold = adaptive_thresh(I,varargin)
% ADPATIVE_THRESH    finds a threshold which separates a given image into
%                    background and signal, assuming a bimodal histrogram
%                    is present, uses the algoritm from p.599-600 in
%                    Digital Image Processing by Gonzalez and Woods
%
%   adaptive_thresh(I) using the pixel values in image 'I', find and return
%   a threshold which should separate the image into background and signal 
%
%   adaptive_thresh(I,W) using the pixel values in image 'I', find and return
%   a threshold which should separate the image into background and signal
%   using a fractional weight 'W' on the upper mean


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upper_mean_weight = 0.5;

if (size(varargin,1) > 0)
    if (isnumeric(varargin{1}))
        if (varargin{1} >= 0 && varargin{1} <= 1)    
            upper_mean_weight = varargin{1};
        else
            error('ERROR: adaptive_thresh - Expected weight value to be between zero and one');
        end
    else
        error('ERROR: adaptive_thresh - Expected second parameter to be a number, the weight on the upper threshold');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = mean(I(:));

threshold_change_thresh = 0.0001;

upper_mean = mean(I(find(I > threshold)));
lower_mean = mean(I(find(I <= threshold)));

new_thresh = upper_mean_weight*upper_mean + (1-upper_mean_weight)*lower_mean;

while (abs(threshold - new_thresh) > threshold_change_thresh)
    threshold = new_thresh;
    
    upper_mean = mean(I(find(I > threshold)));
    lower_mean = mean(I(find(I <= threshold)));

    new_thresh = upper_mean_weight*upper_mean + (1-upper_mean_weight)*lower_mean;
end

end
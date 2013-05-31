Modifying Registration Code

The registration code used at the beginning of the processing pipeline uses a
builtin MATLAB function, but the default matlab implementation doesn't return
the transformation matrix. I need the transformation matrix, so that it can be
applied to both the gel and actin images. So I've modified the code provided
with matlab to return that matrix. The changes are simple, but first you have to
find your copy of imregister.m. I use linux for development, so my copy was
located in /usr/local/MATLAB/R2012a/toolbox/images/images.

After finding the imregister command, you need to copy it into:

src/find_cell_features/matlab_scripts/imregister_trans_out.m

and also copy the folder named 'private' from the MATLAB source directory you
found above.

Now change the first line of imregister_trans_out.m to:

function [moving_reg,t] = imregister_trans_out(varargin)

The only change was to the output variables, adding t, the transformation
matrix.

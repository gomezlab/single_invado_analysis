function [B,F,T] = otsuThresholding( I )

% first get a histogram and normalize it to get the PDF
% use 256 levels

[pdfI,x] = hist( double(I(:)), 256 );
pdfI = pdfI/sum(pdfI);

% we start with all foreground and move to all background

varB = zeros( 256, 1 );
nb = zeros( 256, 1 );
nf = zeros( 256, 1 );
nf(1) = 1;

mub = zeros( 256, 1 );
muf = zeros( 256, 1 );
muf(1) = mean(I(:));

for iI=2:256
  
  nt = pdfI(iI-1);
  
  nb(iI) = nb(iI-1) + nt;
  nf(iI) = nf(iI-1) - nt;
  T = x(iI-1);
  
  mub(iI) = (mub(iI-1)*nb(iI-1)+nt*T)/(nb(iI));
  muf(iI) = (muf(iI-1)*nf(iI-1)-nt*T)/(nf(iI));
  
  varB(iI) = nb(iI)*nf(iI)*(muf(iI)-mub(iI))^2;
  
end

% and now let's find the maximum

[y,indx] = max( varB );

T = x( indx-1 );
B = I;
F = I;

indxB = find( I<T );
indxF = find( I>=T );

B(indxF) = 0;
B(indxB) = 1;

F(indxB) = 0;
F(indxF) = 1;
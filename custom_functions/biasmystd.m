function y = biasmystd(x,blk)
% Hopefully faster, more memory efficient std - Only works for 2D arrays,
% works *across* rows (not down columns)
n = size(x,2);
nRow=size(x,1);
nBlk=nRow/blk;
y = zeros(nRow,1);
I0 = 1:blk;
for i=1:nBlk
  I = I0+(i-1)*blk;
  xbar = sum(x(I,:), 2) ./ n;
  xc = bsxfun(@minus, x(I,:), xbar);
  y(I) = sqrt(sum(xc.^2, 2) ./ n);
end

return

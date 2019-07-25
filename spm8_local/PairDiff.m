function PairDiff(Imgs,BaseNm)
% PairDiff(Imgs,BaseNm)
% Imgs   - Matrix of (even-number of) filenames
% BaseNm - Basename for difference images
%
% Create pairwise difference for a set of images (img2-img1, img4-img3, etc),
% named according to BaseNm (BaseNm_0001, BaseNm_0002, BaseNm_0003, etc).
%______________________________________________________________________
% $Id: PairDiff.m,v 1.2 2012/02/16 21:36:41 nichols Exp $


if nargin<1 | isempty(Imgs)
  Imgs = spm_select(Inf,'image','Select n*2 images');
end
if nargin<2 | isempty(BaseNm)
  BaseNm = spm_input('Enter difference basename','+1','s','Diff');
end
V = spm_vol(Imgs);
n = length(V);
if rem(n,2)~=0
  error('Must specify an even number of images')
end


V1 = V(1);
V1  = rmfield(V1,{'fname','descrip','n','private'});
V1.dt = [spm_type('float32') spm_platform('bigend')];

for i=1:n/2
  V1.fname = sprintf('%s_%04d.img',BaseNm,i);
  V1.descrip = sprintf('Difference %d',i);
  Vo(i) = spm_create_vol(V1);
end


%
% Do the differencing
%

fprintf('Creating differenced data ')

for i=1:2:n
  
  img1 = spm_read_vols(V(i));
  img2 = spm_read_vols(V(i+1));
%   if i==1,
%     refimg = (img1+img2)/2;
%   end

  img = img2-img1; % + refimg;

  spm_write_vol(Vo((i+1)/2),img);

end;


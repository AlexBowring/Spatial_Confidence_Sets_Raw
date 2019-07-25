function [fwhm,VRpv] = spm_est_smoothness(varargin)
% Estimation of smoothness based on [residual] images
% FORMAT [fwhm,VRpv] = spm_est_smoothness(VResI,[VM,df]);
%
% P     - filenames of [standardized residual] images
% PM    - filename of mask image
%
% fwhm  - estimated FWHM in all image directions
% VRpv  - handle of Resels per Voxel image
%_______________________________________________________________________
%  
% spm_est_smoothness returns a spatial smoothness estimator based on the
% variances of the normalized spatial derivatives as described in K.
% Worsley, (1996). Inputs are a mask image and a number of [residual]
% images. Output is a global estimate of the smoothness expressed as the
% FWHM of an equivalent Gaussian point spread function. An estimate of
% resels per voxels (see spm_spm) is written as an image file ('RPV.img')
% to the current directory.
%
% The mask image specifies voxels, used in smoothness estimation, by
% assigning them non-zero values. The dimensions, voxel sizes, orientation 
% of all images must be the same. The dimensions of the images can be of
% dimensions 0, 1, 2 and 3.
% 
% Note that 1-dim images (lines) must exist in the 1st dimension and
% 2-dim images (slices) in the first two dimensions. The estimated fwhm
% for any non-existing dimension is infinity.
%
% 
% Ref:
% 
% K. Worsley (1996). An unbiased estimator for the roughness of a
% multivariate Gaussian random field. Technical Report, Department of
% Mathematics and Statistics, McGill University
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_est_smoothness.m 1143 2008-02-07 19:33:33Z spm $


% assign input argumants
%-----------------------------------------------------------------------
if nargin > 0, V  = varargin{1}; end
if nargin > 1, VM = varargin{2}; end
if nargin > 2, df = varargin{3}; else df = NaN; end
if nargin > 3
    spm('alert!', 'spm_est_smoothness: Wrong number of arguments');
    return;
end
if nargin < 1
    V   = spm_select(inf, '^ResI.*\.img$', 'Select residual images');
end
if nargin < 2
    VM  = spm_select(1, 'mask.img', 'Select mask image');
end

if isnan(df)
  % If no DF set, just assume full
  df=length(V);
end


% intialise
%-----------------------------------------------------------------------
if ~isstruct(V)
    V     = spm_vol(V);
end
if ~isstruct(VM)
    VM    = spm_vol(VM);
end

%-Intialise RESELS per voxel image
%-----------------------------------------------------------------------
VRpv  = struct('fname','RPV.img',...
            'dim',      VM.dim(1:3),...
            'dt',       [spm_type('float64') spm_platform('bigend')],...
            'mat',      VM.mat,...
            'pinfo',    [1 0 0]',...
            'descrip',  'spm_spm: resels per voxel');
VRpv  = spm_create_vol(VRpv);


% dimensionality of image
%-----------------------------------------------------------------------
D     = 3 - sum(VM.dim(1:3) == 1);
if D == 0
    fwhm = [Inf Inf Inf];
    return
elseif D>3
  error('Sorry... 4D random fields aren''t supported')
end

% find voxels within mask
%-----------------------------------------------------------------------
[x,y] = ndgrid(1:VM.dim(1), 1:VM.dim(2));
I     = []; Ix = []; Iy = []; Iz = [];
for i = 1:VM.dim(3)
    z  = i*ones(size(x));
    d  = spm_sample_vol(VM, x, y, z, 0);
    I  = find(d);
    Ix = [Ix; x(I)];
    Iy = [Iy; y(I)];
    Iz = [Iz; z(I)];
end

% compute variance of normalized derivatives in all directions
%-----------------------------------------------------------------------
str   = 'Spatial non-sphericity (over scans)';
fprintf('%-40s: %30s',str,'...estimating derivatives')      %-#
spm_progress_bar('Init',100,'smoothness estimation','');

v     = zeros(size(Ix,1),D);
vO    = zeros(size(Ix,1),D*(D-1)/2);
for i = 1:length(V) % for all residual images
    
    [d, dx, dy, dz] = spm_sample_vol(V(i), Ix, Iy, Iz, 1);
    
    if D >= 1. v(:, 1) = v(:, 1) + dx.^2;  end
    if D >= 2. v(:, 2) = v(:, 2) + dy.^2;  end
    if D >= 3, v(:, 3) = v(:, 3) + dz.^2;  end
    if D == 2
      vO(:,1) = vO(:,1) + dx.*dy;
    elseif D == 3
      % Reverse order helps with determinant calculation below
      vO(:,3) = vO(:,3) + dx.*dy;
      vO(:,2) = vO(:,2) + dx.*dz;
      vO(:,1) = vO(:,1) + dy.*dz;
    end
      
    spm_progress_bar('Set',100*i/length(V));

end
spm_progress_bar('Clear')

% normalise derivatives
%-----------------------------------------------------------------------
v  = v/df;
vO = vO/df;

% eliminate zero variance voxels
%-----------------------------------------------------------------------
I      = find(any(isnan(v')));
v(I,:) = []; Ix(I) = []; Iy(I) = []; Iz(I) = [];


if D == 1
  DL  = v(:,1);
else D == 2
  DL  = prod(v,2) - vO(:,1).^2;
else D == 3
  DL  = prod(v,2) + 2*prod(vO,2) - sum(v.*vO.*vO,2);
end
% Scalar RPV
RPV = (4*log(2))^(-D/2) * sqrt(DL)/sf_DLbias(1/2,D,df);
% D-dim FWHM
FWHM = %%%%% Don't understand what I should do for FWHM... seems like a
       %%%%% nuisance, but infact its key for stationary RFT results.
       %%%%% Buut, does this mean that spm_resel_vol.c is wrong?

% Determinant Lambda to the Half power
DLH 


% resels per voxel (resel) 
% resels = resels/voxel = 1/prod(FWHM)
% FWHM   = sqrt(4.ln2/|dv/dx|))
% fwhm   = 1/FWHM
%-----------------------------------------------------------------------
fprintf('\r%-40s: %30s\n',str,'...writing resels/voxel image')  %-#

fwhm   = sqrt(v./(4*log(2)));
resel  = prod(fwhm,2);
for  i = 1:VM.dim(3)
    d  = NaN*ones(VM.dim(1:2));
    I  = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(VM.dim(1:2), Ix(I), Iy(I))) = resel(I);
    end
    VRpv = spm_write_plane(VRpv, d, i);
end

% global equivalent FWHM {prod(1/FWHM) = (unbiased) RESEL estimator}
%-----------------------------------------------------------------------
fwhm   = mean(fwhm );
RESEL  = mean(resel);
fwhm   = fwhm*((RESEL/prod(fwhm)).^(1/D));
FWHM   = 1./fwhm;
fwhm   = [Inf Inf Inf];
fwhm(1:length(FWHM)) = FWHM;




function c = sf_DLbias(alpha,D,df)
% Correction factor c such that 
%   |Var((d/dx)(e))|^alpha / c
% is unbiased for
%   |Lambda|^alpha
% 
% See Appendix A, S. Hayasaka et al., "Nonstationary cluster-size
% inference with random field and permutation methods", NeuroImage 22
% (2004) 676-687. 

a = 1:D;

c = gamma(df/2-alpha*D)/gamma(df/2) * prod(gamma((df-a)/2+alpha)./gamma((df-a)/2));

return

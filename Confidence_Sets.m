function [] = Confidence_Sets(4D_COPES, GROUP_MASK_IMAGE, THRESH, OUT)

% Computes Confidence Sets and point estimate (yellow) image for subject-level effect estimates maps.  
% Inputs:
%   4D_COPES:         A 4D volume containing all individuals (3D) effect estimates (%BOLD) images. 
%   GROUP_MASK_IMAGE: A 3D image of the group mask
%   THRESH:           The threshold c, in raw change units, mu. 
%   OUT:              Output directory where all output images are saved
%
% Outputs:
%   - Lower_CS.nii        Lower Confidence Set image. All voxels in this image we can assert have an effect size less than c.
%                         In the manuscript, this image was displayed in blue. 
%   - Upper_CS.nii        Upper Confidence Set image. All voxels in this image we can assert have an effect size greater than c.
%                         In the manuscript, this image was displayed in red.
%   - Point_Estimate.nii  Point Estimate Set image. The best guess from the data whre there is an effect size greater than c.
%                         In the manuscript, this image was displayed in yellow.
%   - mean.nii            The sample mean of the effect size images.
%   - sd.nii              The standard deviation of the effect size images. 

tic
[x,y,z] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(x.^2 + y.^2 + z.^2) <=1);

thr   = THRESH;  % In raw change units, mu
nBoot = 5000;

VY=spm_vol(4D_COPES);      % This is the file "handle" for all input
                                 % images  - Ignore gzip warning
VM=spm_vol(GROUP_MASK_IMAGE); % This is the handle for the mask

Mask=spm_read_vols(VM)>0;

nSubj=length(VY);
dim=VY(1).dim;

datamat         = zeros([dim nSubj]);
supG            = zeros([nBoot 1]);
tau             = 1/sqrt(nSubj);

% Load up each of the 3D images into a massive 4D image
for i=1:nSubj
  % Take care to mask each image as its read
  datamat(:,:,:,i) = spm_read_vols(VY(i)).*Mask;
end

observed_mean = mean(datamat,4);
observed_std = std(datamat,0,4); 

observed_AC = observed_mean >= thr;

% Making the interpolated boundary edges
% Horizontal edges
observed_horz = observed_AC(:,2:end,:) | observed_AC(:,1:end-1,:);
% Compute the left shifted horizontal edges
observed_lshift               = observed_AC; % initialize
observed_lshift(:,1:end-1,:)  = observed_horz;
observed_lshift               = observed_lshift & ~observed_AC;
%%% Compute the right shifted horizontal edges
observed_rshift               = observed_AC; % initialize
observed_rshift(:,2:end,:)    = observed_horz;
observed_rshift               = observed_rshift & ~observed_AC;
% Vertical edges
observed_vert = observed_AC(1:end-1,:,:) | observed_AC(2:end,:,:);
%%% Compute the up shifted horizontal edges
observed_ushift               = observed_AC;
observed_ushift(1:end-1,:,:)  = observed_vert;
observed_ushift               = observed_ushift & ~observed_AC;
%%% Compute the down shifted vertical edges
observed_dshift               = observed_AC;
observed_dshift(2:end,:,:)    = observed_vert;
observed_dshift               = observed_dshift & ~observed_AC;
% Depth edges
observed_depth                = observed_AC(:,:,1:end-1) | observed_AC(:,:,2:end);
%%% Compute the back shifted depth edges
observed_bshift               = observed_AC;
observed_bshift(:,:,1:end-1)  = observed_depth;
observed_bshift               = observed_bshift & ~observed_AC;
%%% Compute the front shifted depth edges
observed_fshift              = observed_AC;
observed_fshift(:,:,2:end)   = observed_depth;
observed_fshift              = observed_fshift & ~observed_AC;

% Computing the weights for the weighted linear boundary
observed_lshift_w1 = abs(observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1],:)) - thr)./abs(observed_mean(observed_lshift) - observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1],:)));
observed_lshift_w2 = abs(observed_mean(observed_lshift) - thr)./abs(observed_mean(observed_lshift) - observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1],:)));

observed_rshift_w1 = abs(observed_mean(observed_rshift(:,[2:dim(2) 1],:)) - thr)./abs(observed_mean(observed_rshift) - observed_mean(observed_rshift(:,[2:dim(2) 1],:)));
observed_rshift_w2 = abs(observed_mean(observed_rshift) - thr)./abs(observed_mean(observed_rshift) - observed_mean(observed_rshift(:,[2:dim(2) 1],:)));

observed_ushift_w1 = abs(observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:,:)) - thr)./abs(observed_mean(observed_ushift) - observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:,:)));
observed_ushift_w2 = abs(observed_mean(observed_ushift) - thr)./abs(observed_mean(observed_ushift) - observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:,:)));

observed_dshift_w1 = abs(observed_mean(observed_dshift([2:dim(1) 1],:,:)) - thr)./abs(observed_mean(observed_dshift) - observed_mean(observed_dshift([2:dim(1) 1],:,:)));
observed_dshift_w2 = abs(observed_mean(observed_dshift) - thr)./abs(observed_mean(observed_dshift) - observed_mean(observed_dshift([2:dim(1) 1],:,:)));

observed_bshift_w1 = abs(observed_mean(observed_bshift(:,:,[dim(3) 1:dim(3)-1])) - thr)./abs(observed_mean(observed_bshift) - observed_mean(observed_bshift(:,:,[dim(3) 1:dim(3)-1])));
observed_bshift_w2 = abs(observed_mean(observed_bshift) - thr)./abs(observed_mean(observed_bshift) - observed_mean(observed_bshift(:,:,[dim(3) 1:dim(3)-1])));

observed_fshift_w1 = abs(observed_mean(observed_fshift(:,:,[2:dim(3) 1])) - thr)./abs(observed_mean(observed_fshift) - observed_mean(observed_fshift(:,:,[2:dim(3) 1])));
observed_fshift_w2 = abs(observed_mean(observed_fshift) - thr)./abs(observed_mean(observed_fshift) - observed_mean(observed_fshift(:,:,[2:dim(3) 1])));


% Residuals - just mean centering 
resid = bsxfun(@minus, reshape(datamat, [prod(dim) nSubj]), reshape(observed_mean,[prod(dim) 1]));
resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*resid;

% Finding the values of the residuals on the interpolated boundary
observed_resid_boundary_values = zeros([size(observed_lshift_w1,1)+size(observed_rshift_w1,1)+size(observed_ushift_w1,1)+size(observed_dshift_w1,1)+size(observed_bshift_w1,1)+size(observed_fshift_w1,1) nSubj]);
tic
for i=1:nSubj
  subject_resid_field = reshape(resid(:,i), [dim 1]);

  observed_lshift_boundary_values = observed_lshift_w1.*subject_resid_field(observed_lshift) + observed_lshift_w2.*subject_resid_field(observed_lshift(:,[dim(2) 1:dim(2)-1],:));
  observed_rshift_boundary_values = observed_rshift_w1.*subject_resid_field(observed_rshift) + observed_rshift_w2.*subject_resid_field(observed_rshift(:,[2:dim(2) 1],:));
  observed_ushift_boundary_values = observed_ushift_w1.*subject_resid_field(observed_ushift) + observed_ushift_w2.*subject_resid_field(observed_ushift([dim(1) 1:dim(1)-1],:,:));
  observed_dshift_boundary_values = observed_dshift_w1.*subject_resid_field(observed_dshift) + observed_dshift_w2.*subject_resid_field(observed_dshift([2:dim(1) 1],:,:));
  observed_bshift_boundary_values = observed_bshift_w1.*subject_resid_field(observed_bshift) + observed_bshift_w2.*subject_resid_field(observed_bshift(:,:,[dim(3) 1:dim(3)-1]));
  observed_fshift_boundary_values = observed_fshift_w1.*subject_resid_field(observed_fshift) + observed_fshift_w2.*subject_resid_field(observed_fshift(:,:,[2:dim(3) 1]));

  observed_resid_boundary_values(:,i) = [observed_lshift_boundary_values; observed_rshift_boundary_values; observed_ushift_boundary_values; observed_dshift_boundary_values; observed_bshift_boundary_values; observed_fshift_boundary_values];
end
toc

% Implementing the Multiplier Boostrap to obtain confidence intervals
tic
for k=1:nBoot 
  % Applying the bootstrap using Rademacher variables (signflips)
  signflips                              = randi(2,[nSubj,1])*2-3;

  % Estimated boundary
  observed_boundary_bootstrap       = observed_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
  observed_boundary_resid_field     = sum(observed_boundary_bootstrap, 2)/sqrt(nSubj); 
  % Re-standardizing by bootstrap standard deviation
  observed_boot_std                 = std(observed_boundary_bootstrap, 0, 2);
  observed_boundary_resid_field     = observed_boundary_resid_field./observed_boot_std; 

  supG(k)                  = max(abs(observed_boundary_resid_field));
  
end
toc

supGa95         = prctile(supG,95);


LowerCon  = observed_mean >= thr - supGa95*tau*observed_std;
MiddleCon = observed_AC;
UpperCon  = observed_mean >= thr + supGa95*tau*observed_std;


% Making the edge image for visualization purposes
observed_AC_dil = imdilate(observed_AC,se);
observed_AC_ero = imerode(observed_AC,se);
observed_AC_edge = (observed_AC_dil - observed_AC)|(observed_AC - observed_AC_ero);

subplot(2,3,1)
imagesc(observed_mean(:,:,40));axis image; colorbar
subplot(2,3,2)
hist(supG,50)
abline('v',supGa95)
subplot(2,3,3)
imagesc(observed_AC_edge(:,:,40));axis image; colorbar
subplot(2,3,4)
imagesc(LowerCon(:,:,40));axis image; colorbar
subplot(2,3,5)
imagesc(MiddleCon(:,:,40));axis image; colorbar
subplot(2,3,6)
imagesc(UpperCon(:,:,40));axis image; colorbar

cd(OUT);
Vout=VY(1); % clone the first image's handle
Vout.fname = 'Lower_CS.nii'; % crucially, change the file name!
Vout.descrip = 'Lower Confidence Set. All voxels in this set we can assert have an effect size less than c.'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,LowerCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'Point_Estimate.nii'; % crucially, change the file name!
Vout.descrip = 'The Point Estimate Set. The best guess from the data of voxels with an effect size greater than c.'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,MiddleCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'Upper_CS.nii'; % crucially, change the file name!
Vout.descrip = 'Upper Confidence Set. All voxels in this set we can assert have an effect size greater than c.'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,UpperCon);

% Saving mean, SD and cohen's d image
Vout=VY(1); % clone the first image's handle
Vout.fname = 'mean.nii'; % crucially, change the file name!
Vout.descrip = 'Sample mean of effect estimate images'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,observed_mean);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'sd.nii'; % crucially, change the file name!
Vout.descrip = 'Standard deviation of effect estimate images'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,observed_std);

toc
end
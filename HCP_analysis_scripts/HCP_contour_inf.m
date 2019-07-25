function [] = HCP_Contour_Inf(String,Out)
tic
cd(String);
[x,y,z] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(x.^2 + y.^2 + z.^2) <=1);

thr   = 10.0;  % In raw change units, mu
nBoot = 5000;

VY=spm_vol('smooth_copes.nii.gz');      % This is the file "handle" for all input
                                 % images  - Ignore gzip warning
VM=spm_vol('group_mask.nii.gz'); % This is the handle for the mask

Mask=spm_read_vols(VM)>0;

nSubj=length(VY);
dim=VY(1).dim;

datamat         = zeros([dim nSubj]);
supG            = zeros([nBoot 1]);
supG_cohen_d    = zeros([nBoot 1]);
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

%% Computing variables of interest for the Cohen's d case
cohen_d = observed_mean./observed_std;
cohen_d_std      = sqrt(1 + cohen_d.^2/2);

cohen_d_AC       = cohen_d >= thr;

% Residuals using the SNR transformation
cohen_d_resid = ( (bsxfun(@minus, datamat, observed_mean))./observed_std - cohen_d/2.*(((bsxfun(@minus, datamat, observed_mean))./observed_std).^2-1) )...
                    ./ cohen_d_std;
% Calculating boundary edges for cohen d
% Making the interpolated boundary edges
% Horizontal edges
cohen_d_horz                  = cohen_d_AC(:,2:end,:) | cohen_d_AC(:,1:end-1,:);
% Compute the left shifted horizontal edges
cohen_d_lshift               = cohen_d_AC; % initialize
cohen_d_lshift(:,1:end-1,:)  = cohen_d_horz;
cohen_d_lshift               = cohen_d_lshift & ~cohen_d_AC;
%%% Compute the right shifted horizontal edges
cohen_d_rshift               = cohen_d_AC; % initialize
cohen_d_rshift(:,2:end,:)    = cohen_d_horz;
cohen_d_rshift               = cohen_d_rshift & ~cohen_d_AC;
% Vertical edges
cohen_d_vert			      = cohen_d_AC(1:end-1,:,:) | cohen_d_AC(2:end,:,:);
%%% Compute the up shifted horizontal edges
cohen_d_ushift               = cohen_d_AC;
cohen_d_ushift(1:end-1,:,:)  = cohen_d_vert;
cohen_d_ushift               = cohen_d_ushift & ~cohen_d_AC;
%%% Compute the down shifted vertical edges
cohen_d_dshift               = cohen_d_AC;
cohen_d_dshift(2:end,:,:)    = cohen_d_vert;
cohen_d_dshift               = cohen_d_dshift & ~cohen_d_AC;
% Depth edges
cohen_d_depth                 = cohen_d_AC(:,:,1:end-1) | cohen_d_AC(:,:,2:end);
%%% Compute the back shifted depth edges
cohen_d_bshift               = cohen_d_AC;
cohen_d_bshift(:,:,1:end-1)  = cohen_d_depth;
cohen_d_bshift               = cohen_d_bshift & ~cohen_d_AC;
%%% Compute the front shifted depth edges
cohen_d_fshift              = cohen_d_AC;
cohen_d_fshift(:,:,2:end)   = cohen_d_depth;
cohen_d_fshift              = cohen_d_fshift & ~cohen_d_AC;

% Computing the weights for the weighted linear boundary
cohen_d_lshift_w1 = abs(cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)) - thr)./abs(cohen_d(cohen_d_lshift) - cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)));
cohen_d_lshift_w2 = abs(cohen_d(cohen_d_lshift) - thr)./abs(cohen_d(cohen_d_lshift) - cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)));

cohen_d_rshift_w1 = abs(cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)) - thr)./abs(cohen_d(cohen_d_rshift) - cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)));
cohen_d_rshift_w2 = abs(cohen_d(cohen_d_rshift) - thr)./abs(cohen_d(cohen_d_rshift) - cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)));

cohen_d_ushift_w1 = abs(cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)) - thr)./abs(cohen_d(cohen_d_ushift) - cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)));
cohen_d_ushift_w2 = abs(cohen_d(cohen_d_ushift) - thr)./abs(cohen_d(cohen_d_ushift) - cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)));

cohen_d_dshift_w1 = abs(cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)) - thr)./abs(cohen_d(cohen_d_dshift) - cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)));
cohen_d_dshift_w2 = abs(cohen_d(cohen_d_dshift) - thr)./abs(cohen_d(cohen_d_dshift) - cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)));

cohen_d_bshift_w1 = abs(cohen_d(cohen_d_bshift(:,:,[dim(3) 1:dim(3)-1])) - thr)./abs(cohen_d(cohen_d_bshift) - cohen_d(cohen_d_bshift(:,:,[dim(3) 1:dim(3)-1])));
cohen_d_bshift_w2 = abs(cohen_d(cohen_d_bshift) - thr)./abs(cohen_d(cohen_d_bshift) - cohen_d(cohen_d_bshift(:,:,[dim(3) 1:dim(3)-1])));

cohen_d_fshift_w1 = abs(cohen_d(cohen_d_fshift(:,:,[2:dim(3) 1])) - thr)./abs(cohen_d(cohen_d_fshift) - cohen_d(cohen_d_fshift(:,:,[2:dim(3) 1])));
cohen_d_fshift_w2 = abs(cohen_d(cohen_d_fshift) - thr)./abs(cohen_d(cohen_d_fshift) - cohen_d(cohen_d_fshift(:,:,[2:dim(3) 1])));

% Finding the values of the residuals on the interpolated boundary
cohen_d_resid_boundary_values = zeros([size(cohen_d_lshift_w1,1)+size(cohen_d_rshift_w1,1)+size(cohen_d_ushift_w1,1)+size(cohen_d_dshift_w1,1)+size(cohen_d_bshift_w1,1)+size(cohen_d_fshift_w1,1) nSubj]);
tic
for i=1:nSubj
  subject_resid_field = cohen_d_resid(:,:,:,i);

  cohen_d_lshift_boundary_values = cohen_d_lshift_w1.*subject_resid_field(cohen_d_lshift) + cohen_d_lshift_w2.*subject_resid_field(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:));
  cohen_d_rshift_boundary_values = cohen_d_rshift_w1.*subject_resid_field(cohen_d_rshift) + cohen_d_rshift_w2.*subject_resid_field(cohen_d_rshift(:,[2:dim(2) 1],:));
  cohen_d_ushift_boundary_values = cohen_d_ushift_w1.*subject_resid_field(cohen_d_ushift) + cohen_d_ushift_w2.*subject_resid_field(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:));
  cohen_d_dshift_boundary_values = cohen_d_dshift_w1.*subject_resid_field(cohen_d_dshift) + cohen_d_dshift_w2.*subject_resid_field(cohen_d_dshift([2:dim(1) 1],:,:));
  cohen_d_bshift_boundary_values = cohen_d_bshift_w1.*subject_resid_field(cohen_d_bshift) + cohen_d_bshift_w2.*subject_resid_field(cohen_d_bshift(:,:,[dim(3) 1:dim(3)-1]));
  cohen_d_fshift_boundary_values = cohen_d_fshift_w1.*subject_resid_field(cohen_d_fshift) + cohen_d_fshift_w2.*subject_resid_field(cohen_d_fshift(:,:,[2:dim(3) 1]));

  cohen_d_resid_boundary_values(:,i) = [cohen_d_lshift_boundary_values; cohen_d_rshift_boundary_values; cohen_d_ushift_boundary_values; cohen_d_dshift_boundary_values; cohen_d_bshift_boundary_values; cohen_d_fshift_boundary_values];
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
  
  % cohens d boundary
  %cohen_d_boundary_bootstrap        = cohen_d_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
  %cohen_d_boundary_resid_field      = sum(cohen_d_boundary_bootstrap, 2)/sqrt(nSubj); 
  %supG_cohen_d(k)                   = max(abs(cohen_d_boundary_resid_field));
  
end
toc
supGa95         = prctile(supG,95);
%supGa95_cohen_d = prctile(supG_cohen_d,95);

LowerCon  = observed_mean >= thr - supGa95*tau*observed_std;
MiddleCon = observed_AC;
UpperCon  = observed_mean >= thr + supGa95*tau*observed_std;

%LowerCon_cohen_d   = cohen_d >= thr - supGa95_cohen_d*tau*cohen_d_std;
%MiddleCon_cohen_d  = cohen_d_AC;
%UpperCon_cohen_d   = cohen_d >= thr + supGa95_cohen_d*tau*cohen_d_std;

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

cd(Out);
Vout=VY(1); % clone the first image's handle
Vout.fname = 'smooth_LowerConfidenceInterval_c1000.nii'; % crucially, change the file name!
Vout.descrip = 'Lower confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,LowerCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'smooth_MiddleConfidenceInterval_c1000.nii'; % crucially, change the file name!
Vout.descrip = 'Middle confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,MiddleCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'smooth_UpperConfidenceInterval_c1000.nii'; % crucially, change the file name!
Vout.descrip = 'Upper confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,UpperCon);

% Saving mean, SD and cohen's d image
Vout=VY(1); % clone the first image's handle
Vout.fname = 'smooth_mean.nii'; % crucially, change the file name!
Vout.descrip = 'Sample mean'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,observed_mean);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'smooth_sd.nii'; % crucially, change the file name!
Vout.descrip = 'Standard deviation'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,observed_std);

% Vout=VY(1); % clone the first image's handle
% Vout.fname = 'cohens_d.nii'; % crucially, change the file name!
% Vout.descrip = 'cohens d (signal-to-noise)'; % Actually, put something more
                                        % informative here

% Vout=spm_write_vol(Vout,cohen_d);

%% Creating Cohen d images 
% Vout=VY(1); % clone the first image's handle
% Vout.fname = 'cohen_LowerConfidenceInterval_c0075.nii'; % crucially, change the file name!
% Vout.descrip = 'Lower confidence interval!'; % Actually, put something more
                                        % informative here

% Vout=spm_write_vol(Vout,LowerCon_cohen_d);

% Vout=VY(1); % clone the first image's handle
% Vout.fname = 'cohen_MiddleConfidenceInterval_c0075.nii'; % crucially, change the file name!
% Vout.descrip = 'Middle confidence interval!'; % Actually, put something more
                                        % informative here

% Vout=spm_write_vol(Vout,MiddleCon_cohen_d);

% Vout=VY(1); % clone the first image's handle
% Vout.fname = 'cohen_UpperConfidenceInterval_c0075.nii'; % crucially, change the file name!
% Vout.descrip = 'Upper confidence interval!'; % Actually, put something more
                                        % informative here

% Vout=spm_write_vol(Vout,UpperCon_cohen_d);

toc
end

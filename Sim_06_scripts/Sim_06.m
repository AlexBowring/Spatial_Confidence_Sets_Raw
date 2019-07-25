function Sim_06(nSubj,SvNm,nRlz)

%
% Dimension: 3D
% Signal: Small Sphere (Signal 1)
% Noise: Heterogeneous
%


%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  SvNm  = 'Normsim';  % Save name
end
if (nargin<3)
  nRlz = 5000;
end  
if exist([SvNm '.mat'], 'file')
  error('Will not overwrite sim result')
end

%------------Define parameters
% SvNm = 'LinearSig';
% nSubj  = 120;
% nRlz = 300;

tau     = 1/sqrt(nSubj);
nBoot   = 5000;
dim     = [100 100 100]; 
mag     = 3;
smo     = 3;
rimFWHM = 2/sqrt(2*log(2)); 				 
thr     = 2;
rad     = 5;

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo)*ones(1,3);  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trunc_z     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(3))};
trnind      = cat(2, trunc_x, trunc_y, trunc_z);

resid  = zeros([prod(dim) nSubj]);

% This stores the vector SupG for each run
% This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
subset_success_vector_raw_80           = zeros(nRlz, 1); 
subset_success_vector_raw_90           = zeros(nRlz, 1);
subset_success_vector_raw_95           = zeros(nRlz, 1);
subset_success_vector_observed_80           = zeros(nRlz, 1); 
subset_success_vector_observed_90           = zeros(nRlz, 1);
subset_success_vector_observed_95           = zeros(nRlz, 1);
subset_success_vector_raw_80_alternate = zeros(nRlz, 1); 
subset_success_vector_raw_90_alternate = zeros(nRlz, 1);
subset_success_vector_raw_95_alternate = zeros(nRlz, 1);
subset_success_vector_observed_80_alternate = zeros(nRlz, 1); 
subset_success_vector_observed_90_alternate = zeros(nRlz, 1);
subset_success_vector_observed_95_alternate = zeros(nRlz, 1);

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

threshold_observed_80_store                  = zeros(nRlz, 1);
threshold_observed_90_store                  = zeros(nRlz, 1);
threshold_observed_95_store                  = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store                   = zeros(nBoot, nRlz);
supG_observed_store              = zeros(nBoot, nRlz);

supG_raw                         = zeros(nBoot,1);
supG_observed                    = zeros(nBoot,1);

% Creating a sphere of signal
Sig = SpheroidSignal(wdim, rad, mag, 0);

% Smoothing the signal
Sigs = zeros(wdim);
ss   = spm_smooth(Sig,Sigs,smo*ones(1,3));
Sigs = Sigs;

% Truncate to avoid edge effects
tSigs          = Sigs(trnind{1}, trnind{2}, trnind{3});
maxtSigs       = max(tSigs(:));
Sig            = (mag/maxtSigs)*tSigs;

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;
middle_contour                = AC;
middle_contour_volume         = sum(middle_contour(:));

% Obtaining the edges for the boundary Sig > 2 using the linear interpolation methods 
  % Making the interpolated boundary edges
  % Horizontal edges
  horz = AC(:,2:end,:) | AC(:,1:end-1,:);
  % Compute the left shifted horizontal edges
  lshift              = AC; % initialize
  lshift(:,1:end-1,:) = horz;
  lshift              = lshift & ~AC;
  %%% Compute the right shifted horizontal edges
  rshift              = AC; % initialize
  rshift(:,2:end,:)   = horz;
  rshift              = rshift & ~AC;
  % Vertical edges
  vert = AC(1:end-1,:,:) | AC(2:end,:,:);
  %%% Compute the right shifted horizontal edges
  ushift              = AC;
  ushift(1:end-1,:,:) = vert;
  ushift              = ushift & ~AC;
  %%% Compute the down shifted vertical edges
  dshift              = AC;
  dshift(2:end,:,:)   = vert;
  dshift              = dshift & ~AC;
  % Depth edges
  depth = AC(:,:,1:end-1) | AC(:,:,2:end);
  %%% Compute the back shifted depth edges
  bshift              = AC;
  bshift(:,:,1:end-1) = depth;
  bshift              = bshift & ~AC;
  %%% Compute the front shifted depth edges
  fshift              = AC;
  fshift(:,:,2:end)   = depth;
  fshift              = fshift & ~AC;
  
% Computing the weights for the weighted linear boundary
lshift_w1 = abs(Sig(lshift(:,[dim(2) 1:dim(2)-1],:)) - thr)./abs(Sig(lshift) - Sig(lshift(:,[dim(2) 1:dim(2)-1],:)));
lshift_w2 = abs(Sig(lshift) - thr)./abs(Sig(lshift) - Sig(lshift(:,[dim(2) 1:dim(2)-1],:)));

rshift_w1 = abs(Sig(rshift(:,[2:dim(2) 1],:)) - thr)./abs(Sig(rshift) - Sig(rshift(:,[2:dim(2) 1],:)));
rshift_w2 = abs(Sig(rshift) - thr)./abs(Sig(rshift) - Sig(rshift(:,[2:dim(2) 1],:)));

ushift_w1 = abs(Sig(ushift([dim(1) 1:dim(1)-1],:,:)) - thr)./abs(Sig(ushift) - Sig(ushift([dim(1) 1:dim(1)-1],:,:)));
ushift_w2 = abs(Sig(ushift) - thr)./abs(Sig(ushift) - Sig(ushift([dim(1) 1:dim(1)-1],:,:)));

dshift_w1 = abs(Sig(dshift([2:dim(1) 1],:,:)) - thr)./abs(Sig(dshift) - Sig(dshift([2:dim(1) 1],:,:)));
dshift_w2 = abs(Sig(dshift) - thr)./abs(Sig(dshift) - Sig(dshift([2:dim(1) 1],:,:)));

bshift_w1 = abs(Sig(bshift(:,:,[dim(3) 1:dim(3)-1])) - thr)./abs(Sig(bshift) - Sig(bshift(:,:,[dim(3) 1:dim(3)-1])));
bshift_w2 = abs(Sig(bshift) - thr)./abs(Sig(bshift) - Sig(bshift(:,:,[dim(3) 1:dim(3)-1])));

fshift_w1 = abs(Sig(fshift(:,:,[2:dim(3) 1])) - thr)./abs(Sig(fshift) - Sig(fshift(:,:,[2:dim(3) 1])));
fshift_w2 = abs(Sig(fshift) - thr)./abs(Sig(fshift) - Sig(fshift(:,:,[2:dim(3) 1])));

for t=1:nRlz
    fprintf('.');
    observed_mean = zeros(dim);
    observed_std  = zeros(dim);
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        raw_noise = randn(wdim); %- Noise that will be added to the signal 

        %
        % smooth noise  
        %
        Noises = zeros(wdim);
        tt     = spm_smooth(raw_noise,Noises,smo*ones(1,3));
        Noises = Noises/sqrt(tt);      
      
        %
        % Truncate to avoid edge effects
        %
        tNoises = Noises(trnind{1},trnind{2},trnind{3});
        non_stationary_var = zeros([1 dim(2) dim(3)]);
        non_stationary_var(1,:,:) = repmat(linspace(sqrt(0.5), sqrt(1.5))', [1 dim(2) 1])';
        non_stationary_var = repmat(non_stationary_var, [dim(1) 1 1]);
        tNoises = tNoises.*non_stationary_var;
        tImgs = Sig + tNoises; % Creates the true image of smoothed signal + smoothed noise
        
        observed_mean = observed_mean + tImgs;
        observed_std  = observed_std + tImgs.^2;
        
        tImgs = reshape(tImgs, [prod(dim), 1]);      
        resid(:,i) = tImgs; 
        
      end %========== Loop i (subjects)
      
      observed_mean = observed_mean/nSubj;

      observed_std = sqrt(observed_std/nSubj - observed_mean.^2);
       
      % Making the three observed boundaries: dilated boundary, eroded
      % boundary, and dilated - eroded boundary.
      observed_AC = observed_mean >= thr;
      observed_AC_volume = sum(observed_AC(:)); 
      
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
      
      % Residuals
      resid = bsxfun(@minus,resid,reshape(observed_mean, [prod(dim) 1]));
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*resid;
      
      resid_boundary_values = zeros([size(lshift_w1,1)+size(rshift_w1,1)+size(ushift_w1,1)+size(dshift_w1,1)+size(bshift_w1,1)+size(fshift_w1,1) nSubj]);
      observed_resid_boundary_values = zeros([size(observed_lshift_w1,1)+size(observed_rshift_w1,1)+size(observed_ushift_w1,1)+size(observed_dshift_w1,1)+size(observed_bshift_w1,1)+size(observed_fshift_w1,1) nSubj]);
      for i=1:nSubj
          subject_resid_field = reshape(resid(:,i), [dim 1]);
          lshift_boundary_values = lshift_w1.*subject_resid_field(lshift) + lshift_w2.*subject_resid_field(lshift(:,[dim(2) 1:dim(2)-1],:));
          rshift_boundary_values = rshift_w1.*subject_resid_field(rshift) + rshift_w2.*subject_resid_field(rshift(:,[2:dim(2) 1],:));
          ushift_boundary_values = ushift_w1.*subject_resid_field(ushift) + ushift_w2.*subject_resid_field(ushift([dim(1) 1:dim(1)-1],:,:));
          dshift_boundary_values = dshift_w1.*subject_resid_field(dshift) + dshift_w2.*subject_resid_field(dshift([2:dim(1) 1],:,:));
          bshift_boundary_values = bshift_w1.*subject_resid_field(bshift) + bshift_w2.*subject_resid_field(bshift(:,:,[dim(3) 1:dim(3)-1]));
          fshift_boundary_balues = fshift_w1.*subject_resid_field(fshift) + fshift_w2.*subject_resid_field(fshift(:,:,[2:dim(3) 1]));

          observed_lshift_boundary_values = observed_lshift_w1.*subject_resid_field(observed_lshift) + observed_lshift_w2.*subject_resid_field(observed_lshift(:,[dim(2) 1:dim(2)-1],:));
          observed_rshift_boundary_values = observed_rshift_w1.*subject_resid_field(observed_rshift) + observed_rshift_w2.*subject_resid_field(observed_rshift(:,[2:dim(2) 1],:));
          observed_ushift_boundary_values = observed_ushift_w1.*subject_resid_field(observed_ushift) + observed_ushift_w2.*subject_resid_field(observed_ushift([dim(1) 1:dim(1)-1],:,:));
          observed_dshift_boundary_values = observed_dshift_w1.*subject_resid_field(observed_dshift) + observed_dshift_w2.*subject_resid_field(observed_dshift([2:dim(1) 1],:,:));
          observed_bshift_boundary_values = observed_bshift_w1.*subject_resid_field(observed_bshift) + observed_bshift_w2.*subject_resid_field(observed_bshift(:,:,[dim(3) 1:dim(3)-1]));
          observed_fshift_boundary_values = observed_fshift_w1.*subject_resid_field(observed_fshift) + observed_fshift_w2.*subject_resid_field(observed_fshift(:,:,[2:dim(3) 1]));
          
          resid_boundary_values(:,i) = [lshift_boundary_values; rshift_boundary_values; ushift_boundary_values; dshift_boundary_values; bshift_boundary_values; fshift_boundary_balues];
          observed_resid_boundary_values(:,i) = [observed_lshift_boundary_values; observed_rshift_boundary_values; observed_ushift_boundary_values; observed_dshift_boundary_values; observed_bshift_boundary_values; observed_fshift_boundary_values];
      end
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;

          % True boundary
          boundary_bootstrap                = resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          boundary_resid_field              = sum(boundary_bootstrap, 2)/sqrt(nSubj);
          % Re-standardizing by bootstrap standard deviation
          boot_std                          = std(boundary_bootstrap, 0, 2);
          boundary_resid_field              = boundary_resid_field./boot_std;
          
          supG_raw(k)                       = max(abs(boundary_resid_field));
          
          % Estimated boundary
          observed_boundary_bootstrap       = observed_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          observed_boundary_resid_field     = sum(observed_boundary_bootstrap, 2)/sqrt(nSubj); 
          % Re-standardizing by bootstrap standard deviation
          observed_boot_std                 = std(observed_boundary_bootstrap, 0, 2);
          observed_boundary_resid_field     = observed_boundary_resid_field./observed_boot_std;
          
          supG_observed(k)                  = max(abs(observed_boundary_resid_field));
      end
    
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
       
    lower_contour_raw_80             = observed_mean >= thr - supGa_raw_80*tau*observed_std;
    upper_contour_raw_80             = observed_mean >= thr + supGa_raw_80*tau*observed_std;
    lower_contour_raw_80_volume_prct = sum(lower_contour_raw_80(:))/middle_contour_volume;
    upper_contour_raw_80_volume_prct = sum(upper_contour_raw_80(:))/middle_contour_volume;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_90             = observed_mean >= thr - supGa_raw_90*tau*observed_std;
    upper_contour_raw_90             = observed_mean >= thr + supGa_raw_90*tau*observed_std;
    lower_contour_raw_90_volume_prct = sum(lower_contour_raw_90(:))/middle_contour_volume;
    upper_contour_raw_90_volume_prct = sum(upper_contour_raw_90(:))/middle_contour_volume;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;    
    
    lower_contour_raw_95             = observed_mean >= thr - supGa_raw_95*tau*observed_std;
    upper_contour_raw_95             = observed_mean >= thr + supGa_raw_95*tau*observed_std;
    lower_contour_raw_95_volume_prct = sum(lower_contour_raw_95(:))/middle_contour_volume;
    upper_contour_raw_95_volume_prct = sum(upper_contour_raw_95(:))/middle_contour_volume;
    mid_on_upper_raw_95              = upper_contour_raw_95.*middle_contour;
    lower_on_mid_raw_95              = middle_contour.*lower_contour_raw_95;
    upper_subset_mid_raw_95          = upper_contour_raw_95 - mid_on_upper_raw_95;
    mid_subset_lower_raw_95          = middle_contour - lower_on_mid_raw_95;
    
    % Observed boundary
    supGa_observed_80                     = prctile(supG_observed, 80);
    supGa_observed_90                     = prctile(supG_observed, 90);
    supGa_observed_95                     = prctile(supG_observed, 95);
       
    lower_contour_observed_80             = observed_mean >= thr - supGa_observed_80*tau*observed_std;
    upper_contour_observed_80             = observed_mean >= thr + supGa_observed_80*tau*observed_std;
    lower_contour_observed_80_volume_prct = sum(lower_contour_observed_80(:))/middle_contour_volume;
    upper_contour_observed_80_volume_prct = sum(upper_contour_observed_80(:))/middle_contour_volume;
    mid_on_upper_observed_80              = upper_contour_observed_80.*middle_contour;
    lower_on_mid_observed_80              = middle_contour.*lower_contour_observed_80;
    upper_subset_mid_observed_80          = upper_contour_observed_80 - mid_on_upper_observed_80;
    mid_subset_lower_observed_80          = middle_contour - lower_on_mid_observed_80;
    
    lower_contour_observed_90             = observed_mean >= thr - supGa_observed_90*tau*observed_std;
    upper_contour_observed_90             = observed_mean >= thr + supGa_observed_90*tau*observed_std;
    lower_contour_observed_90_volume_prct = sum(lower_contour_observed_90(:))/middle_contour_volume;
    upper_contour_observed_90_volume_prct = sum(upper_contour_observed_90(:))/middle_contour_volume;
    mid_on_upper_observed_90              = upper_contour_observed_90.*middle_contour;
    lower_on_mid_observed_90              = middle_contour.*lower_contour_observed_90;
    upper_subset_mid_observed_90          = upper_contour_observed_90 - mid_on_upper_observed_90;
    mid_subset_lower_observed_90          = middle_contour - lower_on_mid_observed_90;    
    
    lower_contour_observed_95             = observed_mean >= thr - supGa_observed_95*tau*observed_std;
    upper_contour_observed_95             = observed_mean >= thr + supGa_observed_95*tau*observed_std;
    lower_contour_observed_95_volume_prct = sum(lower_contour_observed_95(:))/middle_contour_volume;
    upper_contour_observed_95_volume_prct = sum(upper_contour_observed_95(:))/middle_contour_volume;
    mid_on_upper_observed_95              = upper_contour_observed_95.*middle_contour;
    lower_on_mid_observed_95              = middle_contour.*lower_contour_observed_95;
    upper_subset_mid_observed_95          = upper_contour_observed_95 - mid_on_upper_observed_95;
    mid_subset_lower_observed_95          = middle_contour - lower_on_mid_observed_95;

    %
    % Storing all variables of interest
    %
    % True boundary variables
    supG_raw_store(:,t)                                    = supG_raw;
    threshold_raw_80_store(t)                              = supGa_raw_80;
    lower_contour_raw_80_volume_prct_store(t)              = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)              = upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                              = supGa_raw_90;
    lower_contour_raw_90_volume_prct_store(t)              = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)              = upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                              = supGa_raw_95;
    lower_contour_raw_95_volume_prct_store(t)              = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)              = upper_contour_raw_95_volume_prct;
    
    % Observed boundary variables
    supG_observed_store(:,t)                                    = supG_observed;
    threshold_observed_80_store(t)                              = supGa_observed_80;
    lower_contour_observed_80_volume_prct_store(t)              = lower_contour_observed_80_volume_prct;
    upper_contour_observed_80_volume_prct_store(t)              = upper_contour_observed_80_volume_prct;
 
    threshold_observed_90_store(t)                              = supGa_observed_90;
    lower_contour_observed_90_volume_prct_store(t)              = lower_contour_observed_90_volume_prct;
    upper_contour_observed_90_volume_prct_store(t)              = upper_contour_observed_90_volume_prct;

    threshold_observed_95_store(t)                              = supGa_observed_95;
    lower_contour_observed_95_volume_prct_store(t)              = lower_contour_observed_95_volume_prct;
    upper_contour_observed_95_volume_prct_store(t)              = upper_contour_observed_95_volume_prct;
    
    lshift_observed_mean_boundary = lshift_w1.*observed_mean(lshift) + lshift_w2.*observed_mean(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_observed_mean_boundary = rshift_w1.*observed_mean(rshift) + rshift_w2.*observed_mean(rshift(:,[2:dim(2) 1],:));
    ushift_observed_mean_boundary = ushift_w1.*observed_mean(ushift) + ushift_w2.*observed_mean(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_observed_mean_boundary = dshift_w1.*observed_mean(dshift) + dshift_w2.*observed_mean(dshift([2:dim(1) 1],:,:));
    bshift_observed_mean_boundary = bshift_w1.*observed_mean(bshift) + bshift_w2.*observed_mean(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_observed_mean_boundary = fshift_w1.*observed_mean(fshift) + fshift_w2.*observed_mean(fshift(:,:,[2:dim(3) 1]));
    
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the true boundary
    lower_condition_80 = thr - supGa_raw_80*tau*observed_std;
    upper_condition_80 = thr + supGa_raw_80*tau*observed_std;
    lower_condition_90 = thr - supGa_raw_90*tau*observed_std;
    upper_condition_90 = thr + supGa_raw_90*tau*observed_std;
    lower_condition_95 = thr - supGa_raw_95*tau*observed_std;
    upper_condition_95 = thr + supGa_raw_95*tau*observed_std;
    
    lshift_lower_condition_80_boundary = lshift_w1.*lower_condition_80(lshift) + lshift_w2.*lower_condition_80(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_80_boundary = rshift_w1.*lower_condition_80(rshift) + rshift_w2.*lower_condition_80(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_80_boundary = ushift_w1.*lower_condition_80(ushift) + ushift_w2.*lower_condition_80(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_80_boundary = dshift_w1.*lower_condition_80(dshift) + dshift_w2.*lower_condition_80(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_80_boundary = bshift_w1.*lower_condition_80(bshift) + bshift_w2.*lower_condition_80(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_80_boundary = fshift_w1.*lower_condition_80(fshift) + fshift_w2.*lower_condition_80(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_80_boundary = lshift_w1.*upper_condition_80(lshift) + lshift_w2.*upper_condition_80(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_80_boundary = rshift_w1.*upper_condition_80(rshift) + rshift_w2.*upper_condition_80(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_80_boundary = ushift_w1.*upper_condition_80(ushift) + ushift_w2.*upper_condition_80(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_80_boundary = dshift_w1.*upper_condition_80(dshift) + dshift_w2.*upper_condition_80(dshift([2:dim(1) 1],:,:));
    bshift_upper_condition_80_boundary = bshift_w1.*upper_condition_80(bshift) + bshift_w2.*upper_condition_80(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_upper_condition_80_boundary = fshift_w1.*upper_condition_80(fshift) + fshift_w2.*upper_condition_80(fshift(:,:,[2:dim(3) 1]));
    
    lshift_lower_condition_90_boundary = lshift_w1.*lower_condition_90(lshift) + lshift_w2.*lower_condition_90(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_90_boundary = rshift_w1.*lower_condition_90(rshift) + rshift_w2.*lower_condition_90(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_90_boundary = ushift_w1.*lower_condition_90(ushift) + ushift_w2.*lower_condition_90(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_90_boundary = dshift_w1.*lower_condition_90(dshift) + dshift_w2.*lower_condition_90(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_90_boundary = bshift_w1.*lower_condition_90(bshift) + bshift_w2.*lower_condition_90(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_90_boundary = fshift_w1.*lower_condition_90(fshift) + fshift_w2.*lower_condition_90(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_90_boundary = lshift_w1.*upper_condition_90(lshift) + lshift_w2.*upper_condition_90(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_90_boundary = rshift_w1.*upper_condition_90(rshift) + rshift_w2.*upper_condition_90(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_90_boundary = ushift_w1.*upper_condition_90(ushift) + ushift_w2.*upper_condition_90(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_90_boundary = dshift_w1.*upper_condition_90(dshift) + dshift_w2.*upper_condition_90(dshift([2:dim(1) 1],:,:));
    bshift_upper_condition_90_boundary = bshift_w1.*upper_condition_90(bshift) + bshift_w2.*upper_condition_90(bshift(:,:,[dim(3) 1:dim(3)-1],:));
    fshift_upper_condition_90_boundary = fshift_w1.*upper_condition_90(fshift) + fshift_w2.*upper_condition_90(fshift(:,:,[2:dim(3) 1]));
    
    lshift_lower_condition_95_boundary = lshift_w1.*lower_condition_95(lshift) + lshift_w2.*lower_condition_95(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_95_boundary = rshift_w1.*lower_condition_95(rshift) + rshift_w2.*lower_condition_95(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_95_boundary = ushift_w1.*lower_condition_95(ushift) + ushift_w2.*lower_condition_95(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_95_boundary = dshift_w1.*lower_condition_95(dshift) + dshift_w2.*lower_condition_95(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_95_boundary = bshift_w1.*lower_condition_95(bshift) + bshift_w2.*lower_condition_95(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_95_boundary = fshift_w1.*lower_condition_95(fshift) + fshift_w2.*lower_condition_95(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_95_boundary = lshift_w1.*upper_condition_95(lshift) + lshift_w2.*upper_condition_95(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_95_boundary = rshift_w1.*upper_condition_95(rshift) + rshift_w2.*upper_condition_95(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_95_boundary = ushift_w1.*upper_condition_95(ushift) + ushift_w2.*upper_condition_95(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_95_boundary = dshift_w1.*upper_condition_95(dshift) + dshift_w2.*upper_condition_95(dshift([2:dim(1) 1],:,:));
    bshift_upper_condition_95_boundary = bshift_w1.*upper_condition_95(bshift) + bshift_w2.*upper_condition_95(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_upper_condition_95_boundary = fshift_w1.*upper_condition_95(fshift) + fshift_w2.*upper_condition_95(fshift(:,:,[2:dim(3) 1]));
    
    lower_condition_80_success = [lshift_observed_mean_boundary < lshift_lower_condition_80_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_80_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_80_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_80_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_80_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_80_boundary];                                 
    upper_condition_80_success = [lshift_observed_mean_boundary >= lshift_upper_condition_80_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_80_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_80_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_80_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_80_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_80_boundary];
                              
    lower_condition_90_success = [lshift_observed_mean_boundary < lshift_lower_condition_90_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_90_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_90_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_90_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_90_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_90_boundary];                              
    upper_condition_90_success = [lshift_observed_mean_boundary >= lshift_upper_condition_90_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_90_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_90_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_90_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_90_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_90_boundary];
                              
    lower_condition_95_success = [lshift_observed_mean_boundary < lshift_lower_condition_95_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_95_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_95_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_95_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_95_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_95_boundary];                              
    upper_condition_95_success = [lshift_observed_mean_boundary >= lshift_upper_condition_95_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_95_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_95_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_95_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_95_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_95_boundary];
                              
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the observed boundary
    lower_condition_80_observed = thr - supGa_observed_80*tau*observed_std;
    upper_condition_80_observed = thr + supGa_observed_80*tau*observed_std;
    lower_condition_90_observed = thr - supGa_observed_90*tau*observed_std;
    upper_condition_90_observed = thr + supGa_observed_90*tau*observed_std;
    lower_condition_95_observed = thr - supGa_observed_95*tau*observed_std;
    upper_condition_95_observed = thr + supGa_observed_95*tau*observed_std;
    
    lshift_lower_condition_80_observed_boundary = lshift_w1.*lower_condition_80_observed(lshift) + lshift_w2.*lower_condition_80_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_80_observed_boundary = rshift_w1.*lower_condition_80_observed(rshift) + rshift_w2.*lower_condition_80_observed(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_80_observed_boundary = ushift_w1.*lower_condition_80_observed(ushift) + ushift_w2.*lower_condition_80_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_80_observed_boundary = dshift_w1.*lower_condition_80_observed(dshift) + dshift_w2.*lower_condition_80_observed(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_80_observed_boundary = bshift_w1.*lower_condition_80_observed(bshift) + bshift_w2.*lower_condition_80_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_80_observed_boundary = fshift_w1.*lower_condition_80_observed(fshift) + fshift_w2.*lower_condition_80_observed(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_80_observed_boundary = lshift_w1.*upper_condition_80_observed(lshift) + lshift_w2.*upper_condition_80_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_80_observed_boundary = rshift_w1.*upper_condition_80_observed(rshift) + rshift_w2.*upper_condition_80_observed(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_80_observed_boundary = ushift_w1.*upper_condition_80_observed(ushift) + ushift_w2.*upper_condition_80_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_80_observed_boundary = dshift_w1.*upper_condition_80_observed(dshift) + dshift_w2.*upper_condition_80_observed(dshift([2:dim(1) 1],:,:)); 
    bshift_upper_condition_80_observed_boundary = bshift_w1.*upper_condition_80_observed(bshift) + bshift_w2.*upper_condition_80_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_upper_condition_80_observed_boundary = fshift_w1.*upper_condition_80_observed(fshift) + fshift_w2.*upper_condition_80_observed(fshift(:,:,[2:dim(3) 1]));
    
    lshift_lower_condition_90_observed_boundary = lshift_w1.*lower_condition_90_observed(lshift) + lshift_w2.*lower_condition_90_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_90_observed_boundary = rshift_w1.*lower_condition_90_observed(rshift) + rshift_w2.*lower_condition_90_observed(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_90_observed_boundary = ushift_w1.*lower_condition_90_observed(ushift) + ushift_w2.*lower_condition_90_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_90_observed_boundary = dshift_w1.*lower_condition_90_observed(dshift) + dshift_w2.*lower_condition_90_observed(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_90_observed_boundary = bshift_w1.*lower_condition_90_observed(bshift) + bshift_w2.*lower_condition_90_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_90_observed_boundary = fshift_w1.*lower_condition_90_observed(fshift) + fshift_w2.*lower_condition_90_observed(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_90_observed_boundary = lshift_w1.*upper_condition_90_observed(lshift) + lshift_w2.*upper_condition_90_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_90_observed_boundary = rshift_w1.*upper_condition_90_observed(rshift) + rshift_w2.*upper_condition_90_observed(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_90_observed_boundary = ushift_w1.*upper_condition_90_observed(ushift) + ushift_w2.*upper_condition_90_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_90_observed_boundary = dshift_w1.*upper_condition_90_observed(dshift) + dshift_w2.*upper_condition_90_observed(dshift([2:dim(1) 1],:,:));
    bshift_upper_condition_90_observed_boundary = bshift_w1.*upper_condition_90_observed(bshift) + bshift_w2.*upper_condition_90_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_upper_condition_90_observed_boundary = fshift_w1.*upper_condition_90_observed(fshift) + fshift_w2.*upper_condition_90_observed(fshift(:,:,[2:dim(3) 1]));
    
    lshift_lower_condition_95_observed_boundary = lshift_w1.*lower_condition_95_observed(lshift) + lshift_w2.*lower_condition_95_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_lower_condition_95_observed_boundary = rshift_w1.*lower_condition_95_observed(rshift) + rshift_w2.*lower_condition_95_observed(rshift(:,[2:dim(2) 1],:));
    ushift_lower_condition_95_observed_boundary = ushift_w1.*lower_condition_95_observed(ushift) + ushift_w2.*lower_condition_95_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_lower_condition_95_observed_boundary = dshift_w1.*lower_condition_95_observed(dshift) + dshift_w2.*lower_condition_95_observed(dshift([2:dim(1) 1],:,:));
    bshift_lower_condition_95_observed_boundary = bshift_w1.*lower_condition_95_observed(bshift) + bshift_w2.*lower_condition_95_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_lower_condition_95_observed_boundary = fshift_w1.*lower_condition_95_observed(fshift) + fshift_w2.*lower_condition_95_observed(fshift(:,:,[2:dim(3) 1]));
    
    lshift_upper_condition_95_observed_boundary = lshift_w1.*upper_condition_95_observed(lshift) + lshift_w2.*upper_condition_95_observed(lshift(:,[dim(2) 1:dim(2)-1],:));
    rshift_upper_condition_95_observed_boundary = rshift_w1.*upper_condition_95_observed(rshift) + rshift_w2.*upper_condition_95_observed(rshift(:,[2:dim(2) 1],:));
    ushift_upper_condition_95_observed_boundary = ushift_w1.*upper_condition_95_observed(ushift) + ushift_w2.*upper_condition_95_observed(ushift([dim(1) 1:dim(1)-1],:,:));
    dshift_upper_condition_95_observed_boundary = dshift_w1.*upper_condition_95_observed(dshift) + dshift_w2.*upper_condition_95_observed(dshift([2:dim(1) 1],:,:));
    bshift_upper_condition_95_observed_boundary = bshift_w1.*upper_condition_95_observed(bshift) + bshift_w2.*upper_condition_95_observed(bshift(:,:,[dim(3) 1:dim(3)-1]));
    fshift_upper_condition_95_observed_boundary = fshift_w1.*upper_condition_95_observed(fshift) + fshift_w2.*upper_condition_95_observed(fshift(:,:,[2:dim(3) 1]));
    
    lower_condition_80_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_80_observed_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_80_observed_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_80_observed_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_80_observed_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_80_observed_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_80_observed_boundary];
    upper_condition_80_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_80_observed_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_80_observed_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_80_observed_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_80_observed_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_80_observed_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_80_observed_boundary];
                              
    lower_condition_90_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_90_observed_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_90_observed_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_90_observed_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_90_observed_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_90_observed_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_90_observed_boundary];
    upper_condition_90_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_90_observed_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_90_observed_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_90_observed_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_90_observed_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_90_observed_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_90_observed_boundary];
                              
    lower_condition_95_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_95_observed_boundary; ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_95_observed_boundary; ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_95_observed_boundary; ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_95_observed_boundary; ...
                                  bshift_observed_mean_boundary < bshift_lower_condition_95_observed_boundary; ...
                                  fshift_observed_mean_boundary < fshift_lower_condition_95_observed_boundary];
    upper_condition_95_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_95_observed_boundary; ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_95_observed_boundary; ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_95_observed_boundary; ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_95_observed_boundary; ...
                                  bshift_observed_mean_boundary >= bshift_upper_condition_95_observed_boundary; ...
                                  fshift_observed_mean_boundary >= fshift_upper_condition_95_observed_boundary];
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the true boundary in mult. bootstrap
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))==0
      subset_success_vector_raw_80(t) = 1;
      fprintf('raw nominal 80 success! \n');
    else 
      subset_success_vector_raw_80(t) = 0; 
      fprintf('raw nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))==0
      subset_success_vector_raw_90(t) = 1;
      fprintf('raw nominal 90 success! \n');
    else 
      subset_success_vector_raw_90(t) = 0; 
      fprintf('raw nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))==0
      subset_success_vector_raw_95(t) = 1;
      fprintf('raw nominal 95 success! \n');
    else 
      subset_success_vector_raw_95(t) = 0; 
      fprintf('raw nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the observed boundary in mult. bootstrap
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))==0
      subset_success_vector_observed_80(t) = 1;
      fprintf('observed nominal 80 success! \n');
    else 
      subset_success_vector_observed_80(t) = 0; 
      fprintf('observed nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))==0
      subset_success_vector_observed_90(t) = 1;
      fprintf('observed nominal 90 success! \n');
    else 
      subset_success_vector_observed_90(t) = 0; 
      fprintf('observed nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))==0
      subset_success_vector_observed_95(t) = 1;
      fprintf('observed nominal 95 success! \n');
    else 
      subset_success_vector_observed_95(t) = 0; 
      fprintf('observed nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the true boundary
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))+sum(upper_condition_80_success)+sum(lower_condition_80_success)==0
      subset_success_vector_raw_80_alternate(t) = 1;
      fprintf('raw nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_80_alternate(t) = 0; 
      fprintf('raw nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))+sum(upper_condition_90_success)+sum(lower_condition_90_success)==0
      subset_success_vector_raw_90_alternate(t) = 1; 
      fprintf('raw nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_90_alternate(t) = 0; 
      fprintf('raw nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))+sum(upper_condition_95_success)+sum(lower_condition_95_success)==0
      subset_success_vector_raw_95_alternate(t) = 1; 
      fprintf('raw nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_95_alternate(t) = 0; 
      fprintf('raw nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))+sum(upper_condition_80_observed_success)+sum(lower_condition_80_observed_success)==0
      subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('observed nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))+sum(upper_condition_90_observed_success)+sum(lower_condition_90_observed_success)==0
      subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('observed nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))+sum(upper_condition_95_observed_success)+sum(lower_condition_95_observed_success)==0
      subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('observed nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('observed nominal 95 alternate true boundary failure! \n');
    end 
                              
end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);

percentage_success_vector_observed_80                    = mean(subset_success_vector_observed_80, 1);
percentage_success_vector_observed_90                    = mean(subset_success_vector_observed_90, 1);
percentage_success_vector_observed_95                    = mean(subset_success_vector_observed_95, 1);

percentage_success_vector_raw_80_alternate               = mean(subset_success_vector_raw_80_alternate, 1);
percentage_success_vector_raw_90_alternate               = mean(subset_success_vector_raw_90_alternate, 1);
percentage_success_vector_raw_95_alternate               = mean(subset_success_vector_raw_95_alternate, 1);

percentage_success_vector_observed_80_alternate          = mean(subset_success_vector_observed_80_alternate, 1);
percentage_success_vector_observed_90_alternate          = mean(subset_success_vector_observed_90_alternate, 1);
percentage_success_vector_observed_95_alternate          = mean(subset_success_vector_observed_95_alternate, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_observed_80_store threshold_observed_90_store threshold_observed_95_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_observed_80 subset_success_vector_observed_90 subset_success_vector_observed_95 subset_success_vector_raw_80_alternate subset_success_vector_raw_90_alternate subset_success_vector_raw_95_alternate subset_success_vector_observed_80_alternate subset_success_vector_observed_90_alternate subset_success_vector_observed_95_alternate '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_observed_80 percentage_success_vector_observed_90 percentage_success_vector_observed_95 percentage_success_vector_raw_80_alternate percentage_success_vector_raw_90_alternate percentage_success_vector_raw_95_alternate percentage_success_vector_observed_80_alternate percentage_success_vector_observed_90_alternate percentage_success_vector_observed_95_alternate '...
      'supG_raw_store supG_observed_store '...
      'middle_contour_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_observed_80_volume_prct_store lower_contour_observed_90_volume_prct_store lower_contour_observed_95_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_observed_80_volume_prct_store upper_contour_observed_90_volume_prct_store upper_contour_observed_95_volume_prct_store'])

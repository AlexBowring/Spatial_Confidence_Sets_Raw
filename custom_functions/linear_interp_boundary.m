function [bdryValues] = linear_interp_boundary(AC, resid_field, nSubj)

% Interpolate values of random field to boundary of an excursion set using
% 4 -connectivity
% Input:
% AC:   random field over a domain in R^dim, thresholded at some value c, and then binarized.
%	   	
% resid_field: The residual field proposed in SSS we wish to bootstrap and find the maximum of along the boundary. Array of size [dim nSubj].          		   			 
% 
% nSubj : number of subjects		
%
% Output:
% bdryValues: a 2D array containing the values of the random fields on the estimated boundary from the threshold. The values
% are interpolations using a 4-connectivity grid and obtained by averaging the points inside and outside the excursion set.

dim = size(AC);
field_dimension = length(dim);

if field_dimension == 2
  % 2D with 4-connectivity	
  % Making the interpolated boundary edges

  % Horizontal edges
  horz = AC(:,2:end) | AC(:,1:end-1);
  % Compute the left shifted horizontal edges
  lshift            = AC; % initialize
  lshift(:,1:end-1) = horz;
  lshift            = lshift & ~AC;
  lshift_points     = sum(lshift(:));
  % Compute the right shifted horizontal edges
  rshift            = AC; % initialize
  rshift(:,2:end)   = horz;
  rshift            = rshift & ~AC;
  rshift_points     = sum(rshift(:));

  % Vertical edges
  vert = AC(1:end-1,:) | AC(2:end,:);
  % Compute the right shifted horizontal edges
  ushift            = AC;
  ushift(1:end-1,:) = vert;
  ushift            = ushift & ~AC;
  ushift_points     = sum(ushift(:));
  % Compute the down shifted vertical edges
  dshift            = AC;
  % Values of random field on down shifted vertical edges
  dshift(2:end,:)   = vert;
  dshift            = dshift & ~AC;
  dshift_points    = sum(dshift(:));

  % Calculating the values of the residual field along each of the boundary edges
  lshift_boundary_values = (resid_field(repmat(lshift, [1, 1, nSubj])) + resid_field(repmat(lshift(:,[dim(2) 1:dim(2)-1]), [1, 1, nSubj])))/2;
  lshift_boundary_values = reshape(lshift_boundary_values, [lshift_points, nSubj]);
  rshift_boundary_values = (resid_field(repmat(rshift, [1, 1, nSubj])) + resid_field(repmat(rshift(:,[2:dim(2) 1]), [1, 1, nSubj])))/2;
  rshift_boundary_values = reshape(rshift_boundary_values, [rshift_points, nSubj]);
  ushift_boundary_values = (resid_field(repmat(ushift, [1, 1, nSubj])) + resid_field(repmat(ushift([dim(1) 1:dim(1)-1],:), [1, 1, nSubj])))/2;
  ushift_boundary_values = reshape(ushift_boundary_values, [ushift_points, nSubj]);
  dshift_boundary_values = (resid_field(repmat(dshift, [1, 1, nSubj])) + resid_field(repmat(dshift([2:dim(1) 1],:), [1, 1, nSubj])))/2;
  dshift_boundary_values = reshape(dshift_boundary_values, [dshift_points, nSubj]);

  bdryValues = [ lshift_boundary_values; rshift_boundary_values; ushift_boundary_values; dshift_boundary_values ];

end 

if field_dimension == 3

    % 3D with 6-connectivity 
	% Making the interpolated boundary edges
	% Horizontal edges
    horz = AC(:,2:end,:) | AC(:,1:end-1,:);
    % Compute the left shifted horizontal edges
    lshift              = AC; % initialize
    lshift(:,1:end-1,:) = horz;
    lshift              = lshift & ~AC;
    lshift_points       = sum(lshift(:));
    % Compute the right shifted horizontal edges
    rshift            	= AC; % initialize
    rshift(:,2:end,:) 	= horz;
    rshift            	= rshift & ~AC;
    rshift_points       = sum(rshift(:));

    % Vertical edges 
    vert = AC(1:end-1,:,:) | AC(2:end,:,:);
    %%% Compute the up shifted vertical edges
    ushift 				= AC;
    ushift(1:end-1,:,:) = vert;
    ushift 				= ushift & ~AC;
    ushift_points       = sum(ushift(:));
    % Compute the down shifted vertical edges
    dshift 				= AC;
    %%% Values of random field on down shifted vertical edges
    dshift(2:end,:,:)   = vert;
    dshift 				= dshift & ~AC;
    dshift_points	 	= sum(dshift(:));

    % depth edges
    depth = AC(:,:,1:end-1) | AC(:,:,2:end);
    % Compute the back shifted depth edges
    bshift 				= AC;
    bshift(:,:,1:end-1) = depth;
    bshift 				= bshift & ~AC;
    bshift_points	 	= sum(bshift(:));
    % Compute the front shifted depth edges
    fshift 				= AC;
    %%% Values of random field on down shifted vertical edges
    fshift(:,:,2:end)   = depth;
    fshift 				= fshift & ~AC;
    fshift_points		= sum(fshift(:));

    % Calculating the values of the residual field along each of the boundary edges
  	lshift_boundary_values = (resid_field(repmat(lshift, [1, 1, 1, nSubj])) + resid_field(repmat(lshift(:,[dim(2) 1:dim(2)-1],:), [1, 1, 1, nSubj])))/2;
  	lshift_boundary_values = reshape(lshift_boundary_values, [lshift_points, nSubj]);    
    rshift_boundary_values = (resid_field(repmat(rshift, [1, 1, 1, nSubj])) + resid_field(repmat(rshift(:,[2:dim(2) 1],:), [1, 1, 1, nSubj])))/2;
    rshift_boundary_values = reshape(rshift_boundary_values, [rshift_points, nSubj]);
    ushift_boundary_values = (resid_field(repmat(ushift, [1, 1, 1, nSubj])) + resid_field(repmat(ushift([dim(1) 1:dim(1)-1],:,:), [1, 1, 1, nSubj])))/2;
    ushift_boundary_values = reshape(ushift_boundary_values, [ushift_points, nSubj]);
    dshift_boundary_values = (resid_field(repmat(dshift, [1, 1, 1, nSubj])) + resid_field(repmat(dshift([2:dim(1) 1],:,:), [1, 1, 1, nSubj])))/2;
    dshift_boundary_values = reshape(dshift_boundary_values, [dshift_points, nSubj]);
    bshift_boundary_values = (resid_field(repmat(bshift, [1, 1, 1, nSubj])) + resid_field(repmat(bshift(:,:,[dim(3) 1:dim(3)-1]), [1, 1, 1, nSubj])))/2;
    bshift_boundary_values = reshape(bshift_boundary_values, [bshift_points, nSubj]);
    fshift_boundary_values = (resid_field(repmat(fshift, [1, 1, 1, nSubj])) + resid_field(repmat(fshift(:,:,[2:dim(3) 1]), [1, 1, 1, nSubj])))/2;
    fshift_boundary_values = reshape(fshift_boundary_values, [fshift_points, nSubj]);

    % concatinated values of field on the linear edges 
    bdryValues = [ lshift_boundary_values; rshift_boundary_values; ushift_boundary_values; dshift_boundary_values; bshift_boundary_values; fshift_boundary_values ]; 

end
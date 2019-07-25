function varargout = CorrClusTh(SPM,u,alpha,guess,msk)
%  Find the corrected cluster size threshold for a given alpha
%  function [k,Pc] =CorrClusTh(SPM,u,alpha,guess)
%  SPM   - SPM data structure
%  u     - Cluster defining threshold
%          If less than zero, u is taken to be uncorrected P-value
%  alpha - FWE-corrected level (defaults to 0.05)
%  guess - Set to NaN to use a Newton-Rhapson search (default)
%          Or provide a explicit list (e.g. 1:1000) of cluster sizes to
%          search over.
%          If guess is a (non-NaN) scalar nothing happens, except the the
%          corrected P-value of guess is printed.
%  msk   - option to provide a mask to restrict the search volume, either
%          as an image volume or as a number (radius of a sphere in mm). 
%
% Finds the corrected cluster size (spatial extent) threshold for a given
% cluster defining threshold u and FWE-corrected level alpha. 
%
%_________________________________________________________________________
% $Id: CorrClusTh.m,v 1.13 2011/12/28 nichols Exp $ Thomas Nichols, Marko Wilke


% settings
epsP = 1e-6;   % Corrected P-value convergence criterion (fraction of alpha)
du   = 1e-6;   % Step-size for Newton-Rhapson
maxi = 100;    % Maximum interations for refined search


% check inputs
if nargin<1 | isempty(SPM)
  load(spm_select(1,'SPM.mat', 'Select SPM.mat to check'));
end

% get data from SPM
  df   = [1 SPM.xX.erdf];                               % degrees of freedom
  STAT = 'T';                                           % assume T-statistics
  n    = 1;                                             % assume no conjunction
  M    = SPM.xVol.M;                                    % get space
  VOX  = sqrt(diag(M(1:3,1:3)'*M(1:3,1:3)))';           % voxel size
  FWHM = SPM.xVol.FWHM;                                 % smoothness {voxels}
  FWHMmm= FWHM.*VOX;                                    % smoothness {mm}
  v2r  = 1/prod(FWHM(~isinf(FWHM)));                    %-voxels to resels


% get initial threshold to use
if nargin<2 | isempty(u)
  u = spm_input(['Cluster defining th. {',STAT,' or p value}'],'+0','r',0.001,1);
end

% if below 1, assume this to be a p-value
if u <= 1; 
  u = spm_u(u,df,STAT); 
end

% get corrected cluster size threshold, assume 0.05, FWE-corrected
if nargin<3 | isempty(alpha)
  alpha = 0.05;
end


% this part by Marko to allow for small-volume correction, using code from spm_VOI and from spm_resels
  if nargin<5 | isempty(msk)

	S    = SPM.xVol.S;
	R    = SPM.xVol.R;
  else

	% user passed a mask, so we should use it
	  if isnumeric(msk)

		% assume that this number represents the radius of a sphere in mm
		  S = msk;
		  S = S./VOX;
		  s = S(:)./FWHM(:);
		  s = s(s > 0);
		  s = prod(s).^(1/3);
		  R = [1 4*s 2*pi*s^2 (4/3)*pi*s^3];
		  S = 4/3*pi*(prod(S).^(1/3))^3;   % volume of the sphere in voxels

	  else

		% what was passed?
		  if ischar(msk)
			V = spm_vol(msk);

		  elseif isstruct(msk)
			V = mask;
		  end;


		% how many voxels in this mask?
		  S = spm_read_vols(V);
		  mVOX   = sqrt(sum(V.mat(1:3,1:3).^2));
		  S = S > 0;
		  S(isnan(S)) = [];
		  S = sum(S(:));


		% also need to recalculate R for this volume
		  FWHM  = FWHM.*(VOX./mVOX);
		  R = spm_resels(FWHM,V,'I');
	  end;
  end;


% show what we got so far
  sf_ShowVolInfo(R,S,Vox,FWHM,FWHMmm)

% get cluster sizes to check, assume brute-force
if nargin<4 | isempty(guess)
  guess = NaN;  % 1:1000;
end


% initialize
epsP = alpha*epsP;
Status = 'OK';

% check options
if length(guess)==1 & ~isnan(guess)
  
  %
  % Dummy case... just report P-value
  %

  k  = guess;
  Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);
  
  Status = 'JustPvalue';

elseif (spm_P(1,1*v2r,u,df,STAT,R,n,S)<alpha)

  %
  % Crazy setting, where 1 voxel cluster is significant
  %

  kr = rad^3;
  k = ceil(kr/v2r);
  Pc  = spm_P(1,k*v2r,u,df,STAT,R,n,S);

  else

	%
	% Brute force!
	%

 	Pc = 1;
	for k = guess

		Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);
		%  fprintf('k=%d Pc=%g\n',k,Pc);
		if Pc <= alpha,  break;  end
	end;
	if (Pc > alpha),  Status = 'OutOfRange';  end
  end


% check, report status
  switch (Status)

	case {'JustPvalue'}

		fprintf(['  For a cluster-defining threshold of %0.4f a cluster size threshold of\n'...
			 '  %d has corrected P-value %g\n\n'],u,k,Pc);

	case {'OK'}

		fprintf(['  For a cluster-defining threshold of %0.4f the level %0.3f corrected\n'...
			 '  cluster size threshold is %d and has size (corrected P-value) %g\n\n'],u,alpha,k,Pc);

	case 'TooRough'

		fprintf(['\n  WARNING: Single voxel cluster is significant!\n\n',...
			 '  For a cluster-defining threshold of %0.4f a k=1 voxel cluster\n'...
			 '  size threshold has size (corrected P-value) %g\n\n'],u,Pc); 

	case 'TooManyIter'

		fprintf(['\n  WARNING: Automated search failed to converge\n' ...
			 '  Try systematic search.\n\n']); 

	case 'OutOfRange'

		fprintf(['\n  WARNING: Within the range of cluster sizes searched (%g...%g)\n',...
			 '  a corrected P-value <= alpha was not found (smallest P: %g)\n\n'],guess(1),guess(end),Pc); 
		fprintf([  '  Try increasing the range or an automatic search.\n\n']); 

	otherwise
		error('Unknown status code');
  end


% generate output
  if nargout==1
	varargout = {k};
  elseif nargout>1
	varargout = {k,Pc};
  end;
  return;
  

function sf_ShowVolInfo(R,S,VOX,FWHM,FWHMmm)

  fprintf('\n  Search Volume:  %7.0f %0.2f x %.02f x %0.2f mm voxels\n',S,VOX);
  fprintf('                    %7.1f cmm, %5.3f L, %5.2f RESELS\n',S*prod(VOX),S*prod(VOX)/100^3,R(end));
  fprintf(['                   %0.2f x %.02f x %0.2f mm FWHM, ','%0.2f x %.02f x %0.2f vox FWHM\n\n'],FWHMmm,FWHM);

return
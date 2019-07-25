function TrSmooth(s,P)
%
% Truncated smoothing
%
%
% $Id: TrSmooth.m,v 1.2 2007/04/06 19:20:53 nichols Exp xiaobih $

if nargin<1, s       = spm_input('smoothing {FWHM in mm}',1); end
if nargin<2, P       = spm_select(Inf,'image','select scans'); end
n     = size(P,1);

for i = 1:n
  Pi = deblank(P(i,:));
  [pth,nm,xt,vr] = fileparts(deblank(Pi));

  % Do standard smoothing, create  's' image
  Pis = fullfile(pth,['s' nm xt vr]);
  spm_smooth(Pi,Pis,s);

  % Create mask, 'm' image
  Pm = fullfile(pth,['m' nm xt vr]);
  spm_imcalc_ui(Pi,Pm,'i1>0',{[],[],spm_type('int16')});

  % Create smoothed mask, 'sm' image
  Pms = fullfile(pth,['sm' nm xt vr]);
  spm_smooth(Pm,Pms,s);

  % Normalize (gain adjust) smoothed image, 'gs'  - tight mask
  Pg = fullfile(pth,['gs' nm xt vr]);
  spm_imcalc_ui({Pis,Pms,Pm},Pg,'i1./i2.*i3',{[],[],spm_type('int16')});

end

 
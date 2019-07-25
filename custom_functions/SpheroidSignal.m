function Img = SpheroidSignal(dim, rad, mag, smo)
%
% function SphereSignal(dim, rad, mag, smo)
% 
% Generate a Ellipsoid of signal with cartesian equasion 
% x^2 + y^2 + z^2 = rad^2
%
% 
% dim = spatial extent of signal
% rad = equatorial radius of Spheroid
% mag = magnitude of signal
% smo = smoothness
%----------------------------------------------
%

dim = dim(:)';
d   = length(dim);

Cent0 = dim/2 + 1/2;
Comm = sprintf('Sphere r=%02d sm=%1d',rad,smo); 
Img = MySmooth(MkRadimg(dim,Cent0)<=rad,smo)*mag;

fprintf('.');

return

function Rmap = MkRadimg(dim,c)
% dim   image dimensions
% c     center
%_________________________________________________________________________


 [x y z] = ndgrid([1:dim(1)],[1:dim(2)],[1:dim(3)]);
 Rmap = sqrt((x-c(1)).^2 + (y-c(2)).^2 + (z-c(3)).^2);


return

function sImg = MySmooth(Img,smo)
% Smooth and normalize to max unity
%_________________________________________________________________________

sImg = zeros(size(Img));

if any(smo)
  spm_smooth(double(Img),sImg,smo);
  sImg = sImg/max(sImg(:));
else
  sImg = Img/max(double(Img(:)));
end

Img(Img(:)<0.05) = 0;

return

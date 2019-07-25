function Img = CircularSignal(Dim, rad, mag, smo)
%
% function CircularSignal(Dim, rad, mag, smo)
% 
% Generate a Circle of signal with cartesian equasion 
% x^2 + y^2 = rad^2
%
% 
% Dim = spatial extent of signal
% rad = radius of circle
% mag = magnitude of signal
% smo = smoothness
%----------------------------------------------
%

Dim = Dim(:)';
d   = length(Dim);

Cent0 = Dim/2 + 1/2;
Comm = sprintf('Circle r=%02d sm=%1d',rad,smo); 
Img = MySmooth(MkRadImg(Dim,Cent0)<=rad,smo)*mag;

fprintf('.');

return





function Rmap = MkRadImg(Dim,c)
% dim   image dimensions
% c     center
%_________________________________________________________________________


 [x,y] = ndgrid([1:Dim(1)],[1:Dim(2)]);
 Rmap = sqrt((x-c(1)).^2 + (y-c(2)).^2);


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

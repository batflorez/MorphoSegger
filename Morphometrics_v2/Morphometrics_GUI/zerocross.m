%Tristan Ursell
%November 2012
%
%Find approximate local zero-crossings in an matrix / image, e.g. finding
%zero crossing of del2 of an image for approximate edge detection
%

function Im_cross=zerocross(Im,varargin)

if ~isempty(varargin)
    sz=abs(varargin{1});
else
    sz=1;
end

%strel1=double(strel('square',sz).getnhood);
strel1=double(strel('disk',sz).getnhood);
Nstrel=sum(strel1(:));

%create signed version of image
Im_sign=sign(Im);
Im_sign(Im_sign==0)=1;

%image normalizer
norm1=conv2(ones(size(Im)),strel1,'same');

%convolve image
Im_conv=conv2(Im_sign,strel1,'same')./norm1*Nstrel;

%create zero crossing image
Im_cross=Im_conv<0;




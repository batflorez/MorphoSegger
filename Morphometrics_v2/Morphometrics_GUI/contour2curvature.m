%Tristan Ursell
%April 2012
%
%kappa=contour2curvature(Xin,Yin);
%
%Xin = input X column vector
%Yin = input Y column vector
%
%kappa =  inverse radius of curvature (pixels) at each point in the wrapped
%contour.
%
%Curvature toward the polygon's center is defined as positive.

function kappa=contour2curvature(Xin,Yin)

%check vector lengths
if length(Xin)~=length(Yin)
    error('Input vectors must be the same length')
end

%check vector lengths
if length(Xin)<5
    error('Input vector lengths must be five or greater.')
end

%correct if not column vectors
if size(Xin,2)>1
    Xin=Xin';
    Yin=Yin';
end

%make sure contour is unwrapped
if and(Xin(1)==Xin(end),Yin(1)==Yin(end))
    Xin=Xin(1:end-1);
    Yin=Yin(1:end-1);
end

%length of inputs
npts=length(Xin);

%generate index list
km=modone((1:npts)-1,npts)';
k0=(1:npts)';
kp=modone((1:npts)+1,npts)';

%calculate tangent vectors
t1=[Xin(k0)-Xin(km),Yin(k0)-Yin(km),zeros(npts,1)];
t2=[Xin(kp)-Xin(k0),Yin(kp)-Yin(k0),zeros(npts,1)];
    
%calculate tangent vector norms
t1norm=sqrt(sum(t1.*t1,2));
t2norm=sqrt(sum(t2.*t2,2));

%calculate cross products
c12=cross(t1/norm(t1),t2/norm(t2));

%calculate the tangent derivative numerator
dtnorm=sqrt(sum((t2-t1).*(t2-t1),2));

%calculate the curvature (1/radius of curvature in pixels)
kappa=sign(c12(:,3)).*dtnorm./(1/2*(t1norm+t2norm));


%MODONE    Modulus after division starting from 1.
function r = modone(x,y)
if y==0
    r = NaN;
else
    n = floor(x/y);
    r = x - n*y;
    r(r==0)=y;
end

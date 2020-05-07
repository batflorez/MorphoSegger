%Tristan Ursell
%Contour Profiler
%July 2015
%

function prof = get_contour_profile(X,Y,im1,dl,res,offset)

%kymograph project length (pixels)
%dl=5;

%profile vector resolution
%res=10;

%offset percentage inside
%offset=0.7;

%initialize rotation matrix
theta0=3/2*pi;
Rotmat=zeros(2,2);
Rotmat(1,1)=cos(theta0);
Rotmat(1,2)=-sin(theta0);
Rotmat(2,1)=sin(theta0);
Rotmat(2,2)=cos(theta0);

%get number of points in contour
npts=length(X);

%****CALCULATING CONTOUR NORMALS AND PROFILES****
prof_cents=zeros(npts,2); %(contour point,x/y)

%point list
pts_list=1:npts;

%point pairs
pt1=modone(pts_list,npts);
pt2=modone(pts_list+1,npts);

%position on contour
prof_cents(:,1)=mean([X(pt1),X(pt2)],2);
prof_cents(:,2)=mean([Y(pt1),Y(pt2)],2);

%calculate tangent vector
n1=[X(pt2)-X(pt1),Y(pt2)-Y(pt1)];
norms=sqrt(n1(:,1).^2 + n1(:,2).^2);
n1=n1./[norms,norms];

%rotate to normal
nrot=(Rotmat*n1')';

%translate and proportion normal vector across the contour
norm1_x=zeros(npts,2);
norm1_y=zeros(npts,2);
norm1_x(:,1)=prof_cents(:,1)-dl*offset*nrot(:,1);
norm1_x(:,2)=prof_cents(:,1)+dl*(1-offset)*nrot(:,1);
norm1_y(:,1)=prof_cents(:,2)-dl*offset*nrot(:,2);
norm1_y(:,2)=prof_cents(:,2)+dl*(1-offset)*nrot(:,2);

%calculate tangents and normals
int_temp=zeros(npts,res*dl);
for k=1:npts
    %calulate projected mean intensity value (kymograph)
    int_temp(k,:)=improfile(im1,norm1_x(k,:),norm1_y(k,:),res*dl,'bilinear')';
end

%perform rank
%pert_keep=0.3;
%int_rank=sort(int_temp,2,'descend');

%save raw kymograph (keep different statistical measures of the
%kyomgraph -- each point on the contour includes a whole line scan
%distribution -- choose your measure)
prof=mean(int_temp,2);


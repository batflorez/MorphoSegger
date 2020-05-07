function res = contour2pill_mesh(coord,stp,tolerance,lng)
% This function performs a medial axis transform to a non-branching axis
% Takes the outline coordinates, step size on the final centerline,
% tolerance to non-centrality, and the length of the ribs. The last
% parameter should be longer than the ribs, but should not be too long to
% intersect the contour again: though most cases of >1 intersection will
% be resolved by the function, some will not. Outputs the coordinates of
% the centerline points.
delta = 1E-10;

%check coordinate inputs
if size(coord,1)<=2
    error('Coordinate inputs too small / wrong dimension for meshing.')
end

%make sure the contour is unwrapped
while true
    if abs(coord(1,1)-coord(end,1))<0.00001 && abs(coord(1,2)-coord(end,2))<0.00001
        coord = coord(1:end-1,:);
    else
        break
    end
end

%get voronoi tesselation
warning('off','MATLAB:delaunay:DuplicateDataPoints')
warning('off','MATLAB:TriRep:PtsNotInTriWarnId')
[vx,vy] = voronoi(coord(:,1),coord(:,2));
warning('on','MATLAB:delaunay:DuplicateDataPoints')
warning('on','MATLAB:TriRep:PtsNotInTriWarnId')

% remove vertices crossing the boundary
vx = reshape(vx,[],1);
vy = reshape(vy,[],1);
q = intxy2(vx,vy,coord([1:end 1],1),coord([1:end 1],2));
vx = reshape(vx(~[q;q]),2,[]);
vy = reshape(vy(~[q;q]),2,[]);

% remove vertices outside
q = logical(inpolygon(vx(1,:),vy(1,:),coord(:,1),coord(:,2))...
    .* inpolygon(vx(2,:),vy(2,:),coord(:,1),coord(:,2)));
vx = reshape(vx([q;q]),2,[]);
vy = reshape(vy([q;q]),2,[]);

% remove isolated points
if isempty(vx), 
    res=0;
    return
end
t = ~((abs(vx(1,:)-vx(2,:))<delta)&(abs(vy(1,:)-vy(2,:))<delta));
vx = reshape(vx([t;t]),2,[]);
vy = reshape(vy([t;t]),2,[]);

% remove branches
vx2=[];
vy2=[];
while true
    for i=1:size(vx,2)
        if((sum(sum((abs(vx-vx(1,i))<delta)&(abs(vy-vy(1,i))<delta)))>1)&&(sum(sum((abs(vx-vx(2,i))<delta)&(abs(vy-vy(2,i))<delta)))>1))
            vx2=[vx2 vx(:,i)];
            vy2=[vy2 vy(:,i)];
        end
    end
    if size(vx,2)-size(vx2,2)<=2,
        vx3 = vx2;
        vy3 = vy2;
        break
    else
        vx = vx2;
        vy = vy2;
        vx2 = [];
        vy2 = [];
    end
end

%{
% check that there are no cycles
vx3 = vx;
vy3 = vy;
vx2 = [];
vy2 = [];
while size(vx,2)>1
    for i=1:size(vx,2)
        if((sum(sum((abs(vx-vx(1,i))<delta)&(abs(vy-vy(1,i))<delta)))>1)&&(sum(sum((abs(vx-vx(2,i))<delta)&(abs(vy-vy(2,i))<delta)))>1))
            vx2=[vx2 vx(:,i)];
            vy2=[vy2 vy(:,i)];
        end
    end
    if size(vx,2)-size(vx2,2)<=1
        res = 0;
        return;
    else
        vx = vx2;
        vy = vy2;
        vx2 = [];
        vy2 = [];
    end
end
vx = vx3;
vy = vy3;
%}

if isempty(vx) || size(vx,1)<2 
    res=0;
    return
end

% sort points
vx2=[];
vy2=[];
for i=1:size(vx,2) % in this cycle find the first point
    if sum(sum(abs(vx-vx(1,i))<delta & abs(vy-vy(1,i))<delta))==1
        vx2=vx(:,i)';
        vy2=vy(:,i)';
        break
    elseif sum(sum(abs(vx-vx(2,i))<delta & abs(vy-vy(2,i))<delta))==1
        vx2=fliplr(vx(:,i)');
        vy2=fliplr(vy(:,i)');
        break
    end
end

k=2;
while true % in this cycle sort all points after the first one
    f1=find(abs(vx(1,:)-vx2(k))<delta & abs(vy(1,:)-vy2(k))<delta & (abs(vx(2,:)-vx2(k-1))>=delta | abs(vy(2,:)-vy2(k-1))>=delta));
    f2=find(abs(vx(2,:)-vx2(k))<delta & abs(vy(2,:)-vy2(k))<delta & (abs(vx(1,:)-vx2(k-1))>=delta | abs(vy(1,:)-vy2(k-1))>=delta));
    if f1>0
        vx2 = [vx2 vx(2,f1)];
        vy2 = [vy2 vy(2,f1)];
    elseif f2>0
        vx2 = [vx2 vx(1,f2)];
        vy2 = [vy2 vy(1,f2)];
    else
        break
    end
    k=k+1;
end

skel=[vx2' vy2'];

if size(vx2,2)<=1
    res = 0;
    return
end

% interpolate skeleton to equal step, extend outside of the cell and smooth
% tolerance=0.001;
d=diff(skel,1,1);
l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
if l(end)>=stp
    skel = [interp1(l,vx2,0:stp:l(end))' interp1(l,vy2,0:stp:l(end))'];
else
    skel = [vx2(1) vy2(1);vx2(end) vy2(end)];
end

if size(skel,1)<=1 % || length(skel)<4
    res = 0;
    return
end

lng0 = l(end);
sz = lng0/stp;
L = size(coord,1);
coord = [coord;coord(1,:)];
% lng = 100;
skel2 = [skel(1,:)*(lng/stp+1) - skel(2,:)*lng/stp; skel;...
    skel(end,:)*(lng/stp+1) - skel(end-1,:)*lng/stp];
d=diff(skel2,1,1);
l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
[l,i] = unique(l);
skel2 = skel2(i,:);

%TSU: 1/7/2015
%Issue with skel2 not reaching the poles and thus no intersection of
%centerline with contour -- this is what leads to error below. Increaasing
%the 'lng' parameter fixes this, but at what cost?

if length(l)<2 || size(skel2,1)<2
    res=0; 
    return
end

% find the intersection of the 1st skel2 segment with the contour, the
% closest to one of the poles (which will be called 'prevpoint')
[~,~,indS,indC]=intxyMulti(skel2(2:-1:1,1),skel2(2:-1:1,2),coord(:,1),coord(:,2));  
[~,prevpoint] = min([min(modL(indC,1)) min(modL(indC,L/2+1))]);
if prevpoint==2
    prevpoint=L/2; 
end

% prevpoint = mod(round(indC(1))-1,L)+1;
skel3=spsmooth(l,skel2',tolerance,0:stp:l(end))';%1:stp:

% recenter and smooth again the skeleton
[pintx,pinty,q]=skel2mesh(skel3);
if length(pintx)<sz-1
    skel3=spsmooth(l,skel2',tolerance/100,0:stp:l(end))';
    [pintx,pinty,q]=skel2mesh(skel3);
end

if ~q || length(pintx)<sz-1  || length(skel3)<4
    res=-1;
    return
end

skel = [mean(pintx,2) mean(pinty,2)];
d=diff(skel,1,1);
l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
[l,i] = unique(l);
skel = skel(i,:);
if length(l)<2 || size(skel,1)<2, res=0; return; end
skel=spsmooth(l,skel',tolerance,-3*stp:stp:l(end)+4*stp)';

% get the mesh
[pintx,pinty,q]=skel2mesh(skel);
if ~q
    res=0;
    return
end

res = [pintx(:,1) pinty(:,1) pintx(:,2) pinty(:,2)];
if (pintx(1,1)-coord(end,1))^2+(pinty(1,1)-coord(end,2))^2 > (pintx(end,1)-coord(end,1))^2+(pinty(end,1)-coord(end,2))^2
    res = flipud(res);
end

if numel(res)==1 || length(res)<=4
    res=0; 
end

if length(res)>1 && (res(1,1)~=res(1,3) || res(end,1)~=res(end,3))
    res=0; 
end

    function out = modL(in,shift)
        out = mod(in-shift,L);
        out = min(out,L-out);
    end

    function [pintx,pinty,q]=skel2mesh(sk)
        % This function finds intersections of ribs with the contour
        % To be used in "model2mesh" function
        
        if isempty(sk), pintx=[]; pinty=[]; q=false; return; end
        % Find the intersection of the skel with the contour closest to prevpoint
        pintx=[];
        pinty=[];
        [intX,intY,indS,indC]=intxyMulti(sk(:,1),sk(:,2),coord(:,1),coord(:,2));
        if isempty(intX) || isempty(indC) || isempty(prevpoint)
            q=false;
            return;
        end
        [prevpoint,ind] = min(modL(indC,prevpoint));
        prevpoint = indC(ind);
        indS=indS(ind);
        if indS>(size(sk,1)+1-indS)
            sk = sk(ceil(indS):-1:1,:);
        else
            sk = sk(floor(indS):end,:);
        end
        % 2. define the first pair of intersections as this point
        % 3. get the list of intersections for the next pair
        % 4. if more than one, take the next in the correct direction
        % 5. if no intersections found in the reqion between points, stop
        % 6. goto 3.
        % Define the lines used to compute intersections
        d=diff(sk,1,1);
        plinesx1 = repmat(sk(1:end-1,1),1,2)+lng/stp*d(:,2)*[0 1];
        plinesy1 = repmat(sk(1:end-1,2),1,2)-lng/stp*d(:,1)*[0 1];
        plinesx2 = repmat(sk(1:end-1,1),1,2)+lng/stp*d(:,2)*[0 -1];
        plinesy2 = repmat(sk(1:end-1,2),1,2)-lng/stp*d(:,1)*[0 -1];
        % Allocate memory for the intersection points
        pintx = zeros(size(sk,1)-1,2);
        pinty = zeros(size(sk,1)-1,2);
        % Define the first pair of intersections as the prevpoint
        pintx(1,:) = [intX(ind) intX(ind)];
        pinty(1,:) = [intY(ind) intY(ind)];
        prevpoint1 = prevpoint;
        prevpoint2 = prevpoint;
        
        % for i=1:size(d,1), plot(plinesx(i,:),plinesy(i,:),'r'); end % Testing
        q=true;
        fg = 1;
        jmax = size(sk,1)-1;
        for j=2:jmax
            % gdisp(['Use 1: ' num2str(size(plinesx1(j,:))) ' ' num2str(size(coord(:,1))) ' j=' num2str(j)]); % Testing
            [pintx1,pinty1,~,indC1]=intxyMulti(plinesx1(j,:),plinesy1(j,:),coord(:,1),coord(:,2),floor(prevpoint1),1);%
            [pintx2,pinty2,~,indC2]=intxyMulti(plinesx2(j,:),plinesy2(j,:),coord(:,1),coord(:,2),ceil(prevpoint2),-1);%
            if (~isempty(pintx1))&&(~isempty(pintx2))
                if pintx1~=pintx2
                    if fg==3
                        break;
                    end
                    fg = 2;
                    [prevpoint1,ind1] = min(modL(indC1,prevpoint1));
                    [prevpoint2,ind2] = min(modL(indC2,prevpoint2));
                    prevpoint1 = indC1(ind1);
                    prevpoint2 = indC2(ind2);
                    pintx(j,:)=[pintx1(ind1) pintx2(ind2)];
                    pinty(j,:)=[pinty1(ind1) pinty2(ind2)];
                else
                    q=false;
                    return
                end
            elseif fg==2
                fg = 3;
            end
        end
        pinty = pinty(pintx(:,1)~=0,:);
        pintx = pintx(pintx(:,1)~=0,:);
        [intX,intY,indS,indC]=intxyMulti(sk(:,1),sk(:,2),coord(:,1),coord(:,2));
        [prevpoint,ind] = max(modL(indC,prevpoint));
        pintx = [pintx;[intX(ind) intX(ind)]];
        pinty = [pinty;[intY(ind) intY(ind)]];
        nonan = ~isnan(pintx(:,1))&~isnan(pinty(:,1))&~isnan(pintx(:,2))&~isnan(pinty(:,2));
        pintx = pintx(nonan,:);
        pinty = pinty(nonan,:);
    end
end

function res=intxy2(ax,ay,bx,by)
% modified intxy, finds the first point along line a where an intersection
% rewritten to accept unsorted sets, 2xN of N segments
% finds all intersections, only gives the row numbers in the first set

try
    
    if isempty(ax)
        res=[];
    else
        res=intxy2C(ax,ay,bx,by);
    end
    
catch
    
    % below is a MATLAB version of the C/C++ routine above, slow but same results
    res = zeros(size(ax,2),1);
    for i=1:size(ax,2)
        %vector to next vertex in a
        u=[ax(2,i)-ax(1,i) ay(2,i)-ay(1,i)];
        %go through each vertex in b
        for j=1:size(bx,2)
            %check for intersections at the vortices
            if (ax(1,i)==bx(1,j) && ay(1,i)==by(1,j)) || (ax(2,i)==bx(1,j) && ay(2,i)==by(1,j))...
                    || (ax(1,i)==bx(2,j) && ay(1,i)==by(2,j)) ||(ax(2,i)==bx(2,j) && ay(2,i)==by(2,j))
                res(i) = 1;
                continue
            end
            %vector from ith vertex in a to jth vertex in b
            v=[bx(2,j)-ax(1,i) by(2,j)-ay(1,i)];
            %vector from ith+1 vertex in a to jth vertex in b
            w=[bx(1,j)-ax(2,i) by(1,j)-ay(2,i)];
            %vector from ith vertex of a to jth-1 vertex of b
            vv=[bx(1,j)-ax(1,i) by(1,j)-ay(1,i)];
            %vector from jth-1 vertex of b to jth vertex of b
            z=[bx(2,j)-bx(1,j) by(2,j)-by(1,j)];
            %cross product of u and v
            cpuv=u(1)*v(2)-u(2)*v(1);
            %cross product of u and vv
            cpuvv=u(1)*vv(2)-u(2)*vv(1);
            %cross product of v and z
            cpvz=v(1)*z(2)-v(2)*z(1);
            %cross product of w and z
            cpwz=w(1)*z(2)-w(2)*z(1);
            % check if there is an intersection
            if cpuv*cpuvv<0 && cpvz*cpwz<0
                res(i) = 1;
            end
        end
    end
end
end

function out = modL(in,shift)
L = 200;
out = mod(in-shift,L);
out = min(out,L-out);
end

function res = spsmooth(x,y,p,xx)
% A simple smoothing cubic spline routine. It takes original points y,
% their parameterization x, tolerance p, and the new parameterization xx.
xi = reshape(x,[],1); % make x vertical
yi = y'; % make y vertival
n = size(xi,1);
ny = size(yi,2);
nn = ones(1,ny);
dx = diff(xi);
drv = diff(yi)./dx(:,nn);
% adx = abs(dx);
% w = max([adx;0],[0;adx])/mean(adx);
w = ones(length(x),1);
if n>2
    idx = 1./dx;
    R = spdiags([dx(2:n-1), 2*(dx(2:n-1)+dx(1:n-2)), dx(1:n-2)], -1:1, n-2, n-2);
    Qt = spdiags([idx(1:n-2), -(idx(2:n-1)+idx(1:n-2)), idx(2:n-1)], 0:2, n-2, n);
    W = spdiags(w,0,n,n);
    Qtw = Qt*spdiags(sqrt(w),0,n,n);
    u = ((6*(1-p))*(Qtw*Qtw.')+p*R)\diff(drv);
    yi = yi - (6*(1-p))*W*diff([zeros(1,ny); diff([zeros(1,ny);u;zeros(1,ny)])./dx(:,nn); zeros(1,ny)]);
    c3 = [zeros(1,ny);p*u;zeros(1,ny)];
    c2 = diff(yi)./dx(:,nn)-dx(:,nn).*(2*c3(1:n-1,:)+c3(2:n,:));
    coefs = reshape([(diff(c3)./dx(:,nn)).',3*c3(1:n-1,:).',c2.',yi(1:n-1,:).'],(n-1)*ny,4);
else % straight line output
    coefs = [drv.' yi(1,:).'];
end
breaks = xi.';
sizec = size(coefs);
k = sizec(end);
l = prod(sizec(1:end-1))/ny;
[mx,nx] = size(xx);
lx = mx*nx;
xs = reshape(xx,1,lx);
[~,index] = histc(xs,[-inf,breaks(2:end-1),inf]);
NaNx = find(index==0);
index = min(index,numel(breaks)-1);
index(NaNx) = 1;
xs = xs-breaks(index);
xs = reshape(repmat(xs,ny,1),1,ny*lx);
index = reshape(repmat(1+ny*index,ny,1)+repmat((-ny:-1).',1,lx), ny*lx, 1 );
v = coefs(index,1).';
for i=2:k
    v = xs.*v + coefs(index,i).';
end
if ~isempty(NaNx) && k==1 && l>1, v = reshape(v,ny,lx); v(:,NaNx) = NaN; end
v = reshape(v,ny*mx,nx);

res = reshape(v,[ny,length(xx)]);
end

function [aintx,ainty,aia,aib]=intxyMulti(ax,ay,bx,by,varargin)
% Finds all intersection points between lines ax,ay and bx,by
% Outputs intersection coordinates aintx,ainty and the position on a and b,
% counted from 1 (first point) to n (last point), real number. All outputs
% are vectors if multiple intersections occurs
% If any other argument provided, only the first intersection point is
% given
% If 2 additional arguments provided, the first one is the starting point,
% the next one (1 or -1) the direction (in the second/"b" variable)

try
    
    if isempty(varargin), jstart = -1; jdir = 1;
    elseif length(varargin)==1, jstart = 1; jdir = 1;
    elseif length(varargin)==2, jstart = varargin{1}; jdir = varargin{2};
    end
    [aintx,ainty,aia,aib]=intxyMultiC(ax,ay,bx,by,jstart,jdir);
    
catch
    
    %Below is a Matlab version of the function (gives the same results)
    jstart = 1;
    jdir = 1;
    if isempty(varargin), multiout = true;
    elseif length(varargin)==1, multiout=false;
    elseif length(varargin)==2, multiout=false; jstart = varargin{1}; jdir = varargin{2};
    end
    jrange = mod((jstart:jdir:length(bx)*jdir+jstart-1)-1,length(bx))+1;
    aintx = [];
    ainty = [];
    aia = [];
    aib = [];
    
    for i=1:length(ax)-1
        %vector to next vertex in a
        u=[ax(i+1)-ax(i) ay(i+1)-ay(i)];
        %normalize u
        magu=sqrt(u(1)^2+u(2)^2);
        if magu==0, continue; end
        u=u/magu;
        %go through each vertex in b
        for j=jrange;%j=1:length(bx)
            %check for intersection
            if ax(i)==bx(j) && ay(i)==by(j)
                aintx = [aintx ax(i)];
                ainty = [ainty ay(i)];
                aia = [aia i];
                aib = [aib j];
                continue
            end
            %vector from ith vertex in a to jth vertex in b
            v=[bx(j)-ax(i) by(j)-ay(i)];
            %vector from ith+1 vertex in a to jth vertex in b
            w=[bx(j)-ax(i+1) by(j)-ay(i+1)];
            %check whether a and b crossed
            if j>1
                %vector from ith vertex of a to jth-1 vertex of b
                vv=[bx(j-1)-ax(i) by(j-1)-ay(i)];
                %vector from jth-1 vertex of b to jth vertex of b
                z=[bx(j)-bx(j-1) by(j)-by(j-1)];
                %cross product of u and v
                cpuv=u(1)*v(2)-u(2)*v(1);
                %cross product of u and vv
                cpuvv=u(1)*vv(2)-u(2)*vv(1);
                %cross product of v and z
                cpvz=v(1)*z(2)-v(2)*z(1);
                %cross product of w and z
                cpwz=w(1)*z(2)-w(2)*z(1);
                if cpuv*cpuvv<=0 && cpvz*cpwz<0
                    %normalize v
                    magv=sqrt(v(1)^2+v(2)^2);
                    if magv==0, continue; end
                    v=v/magv;
                    %normalize z
                    magz=sqrt(z(1)^2+z(2)^2);
                    if magz==0; continue; end
                    z=z/magz;
                    %ith segment of a crosses jth segment of b
                    %cpa=0;
                    %range from a(i) to intersection (Law of Sines)
                    r=magv*sin(acos(dot(v,z)))/sin(acos(dot(-u,z)));
                    if cpuv*cpuvv<=0 && r==1, continue; end
                    %record index along a to include portion between ith and
                    %ith+1 vertex where cpa occurs
                    ia=i+r/magu;
                    %project u along x and y to find point of itersection
                    intx=ax(i)+r/magu*(ax(i+1)-ax(i));
                    inty=ay(i)+r/magu*(ay(i+1)-ay(i));
                    %range from b(j) to intersection
                    r=sqrt((bx(j)-intx)^2+(by(j)-inty)^2);
                    if cpuv*cpuvv<=0 && r==1, continue; end
                    %record index along b to include portion between jth-1 and
                    %jth vertex where intersection occurs
                    ib=j-r/magz;
                    %assign to output variables
                    aintx = [aintx intx];
                    ainty = [ainty inty];
                    aia = [aia ia];
                    aib = [aib ib];
                    if ~multiout, return; end
                end
            end
        end
    end
end
end


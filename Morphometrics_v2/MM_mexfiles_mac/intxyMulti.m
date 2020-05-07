function [aintx,ainty,aia,aib]=intxyMulti(ax,ay,bx,by,varargin)
% Finds all intersection points between lines ax,ay and bx,by
% Outputs intersection coordinates aintx,ainty and the position on a and b,
% counted from 1 (first point) to n (last point), real number. All outputs
% are vectors if multiple intersections occurs
% If any other argument provided, only the first intersection point is
% given
% If 2 additional arguments provided, the first one is the starting point,
% the next one (1 or -1) the direction (in the second/"b" variable)

if isempty(varargin), jstart = -1; jdir = 1;
elseif length(varargin)==1, jstart = 1; jdir = 1;
elseif length(varargin)==2, jstart = varargin{1}; jdir = varargin{2};
end
[aintx,ainty,aia,aib]=intxyMultiC(ax,ay,bx,by,jstart,jdir);

% Below is a Matlab version of the function (gives the same results)

% jstart = 1; 
% jdir = 1;
% if isempty(varargin), multiout = true; 
% elseif length(varargin)==1, multiout=false;
% elseif length(varargin)==2, multiout=false; jstart = varargin{1}; jdir = varargin{2};
% end
% jrange = mod((jstart:jdir:length(bx)*jdir+jstart-1)-1,length(bx))+1;
% aintx = [];
% ainty = [];
% aia = [];
% aib = [];
% 
% for i=1:length(ax)-1
%     %vector to next vertex in a
%     u=[ax(i+1)-ax(i) ay(i+1)-ay(i)];
%     %normalize u
%     magu=sqrt(u(1)^2+u(2)^2);
%     if magu==0, continue; end
%     u=u/magu; 
%     %go through each vertex in b
%     for j=jrange;%j=1:length(bx)
%         %check for intersection
%         if ax(i)==bx(j) && ay(i)==by(j)
%             aintx = [aintx ax(i)];
%             ainty = [ainty ay(i)];
%             aia = [aia i];
%             aib = [aib j];
%             continue
%         end
%         %vector from ith vertex in a to jth vertex in b
%         v=[bx(j)-ax(i) by(j)-ay(i)];
%         %vector from ith+1 vertex in a to jth vertex in b
%         w=[bx(j)-ax(i+1) by(j)-ay(i+1)];
%         %check whether a and b crossed
%         if j>1
%             %vector from ith vertex of a to jth-1 vertex of b
%             vv=[bx(j-1)-ax(i) by(j-1)-ay(i)];
%             %vector from jth-1 vertex of b to jth vertex of b
%             z=[bx(j)-bx(j-1) by(j)-by(j-1)];
%             %cross product of u and v
%             cpuv=u(1)*v(2)-u(2)*v(1);
%             %cross product of u and vv
%             cpuvv=u(1)*vv(2)-u(2)*vv(1);
%             %cross product of v and z
%             cpvz=v(1)*z(2)-v(2)*z(1);
%             %cross product of w and z
%             cpwz=w(1)*z(2)-w(2)*z(1);   
%             if cpuv*cpuvv<=0 && cpvz*cpwz<0
%                 %normalize v
%                 magv=sqrt(v(1)^2+v(2)^2);
%                 if magv==0, continue; end
%                 v=v/magv;
%                 %normalize z
%                 magz=sqrt(z(1)^2+z(2)^2);
%                 if magz==0; continue; end
%                 z=z/magz;  
%                 %ith segment of a crosses jth segment of b
%                 %cpa=0;
%                 %range from a(i) to intersection (Law of Sines)
%                 r=magv*sin(acos(dot(v,z)))/sin(acos(dot(-u,z)));
%                 if cpuv*cpuvv<=0 && r==1, continue; end
%                 %record index along a to include portion between ith and
%                 %ith+1 vertex where cpa occurs
%                 ia=i+r/magu;
%                 %project u along x and y to find point of itersection
%                 intx=ax(i)+r/magu*(ax(i+1)-ax(i));
%                 inty=ay(i)+r/magu*(ay(i+1)-ay(i));
%                 %range from b(j) to intersection 
%                 r=sqrt((bx(j)-intx)^2+(by(j)-inty)^2);
%                 if cpuv*cpuvv<=0 && r==1, continue; end
%                 %record index along b to include portion between jth-1 and
%                 %jth vertex where intersection occurs
%                 ib=j-r/magz;
%                 %assign to output variables
%                 aintx = [aintx intx];
%                 ainty = [ainty inty];
%                 aia = [aia ia];
%                 aib = [aib ib];
%                 if ~multiout, return; end
%             end
%         end
%     end
% end
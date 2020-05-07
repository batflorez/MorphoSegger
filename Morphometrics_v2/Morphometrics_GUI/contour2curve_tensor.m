%Tristan Ursell
%April 2012
%
% [kappa1,kappa2] = contour2curve_tensor(Xcont,Ycont,Xcent,Ycent);
%
%Xcont = contour X column vector
%Ycont = contour Y column vector
%
%Xcent = centerline X column vector
%Ycent = centerline Y column vector
%
%kappa1 & kappa2 =  inverse principle radii of curvature (pixels) at each
%point on the wrapped contour.
%
%Curvature toward the polygon's center is defined as positive.

function [kappa1,kappa2]=contour2curve_tensor(Xcont,Ycont,Xcent,Ycent,maxL)

%check vector lengths
if length(Xcont)~=length(Ycont)
    error('Contour input vectors must be the same length')
end

%check vector lengths
if length(Xcont)<6
    error('Input vector lengths must be six or greater.')
end

%check vector lengths
if length(Xcent)~=length(Ycent)
    error('Contour input vectors must be the same length')
end

%check vector lengths
if length(Xcent)<4
    error('Centerline input vector lengths must be four or greater.')
end

%correct if not column vectors
if size(Xcont,2)>1
    Xcont=Xcont';
    Ycont=Ycont';
end

%correct if not column vectors
if size(Xcent,2)>1
    Xcent=Xcent';
    Ycent=Ycent';
end

%make sure contour is unwrapped
if and(Xcont(1)==Xcont(end),Ycont(1)==Ycont(end))
    Xcont=Xcont(1:end-1);
    Ycont=Ycont(1:end-1);
end

%make sure contour is unwrapped
if and(Xcent(1)==Xcent(end),Ycent(1)==Ycent(end))
    Xcent=Xcent(1:end-1);
    Ycent=Ycent(1:end-1);
end

%calculate first principal curvature
kappa1_other=contour2curvature(Xcont,Ycont);

%length of inputs
npts=length(Xcont);

kappa1=zeros(npts,1);
kappa2=zeros(npts,1);

%generate index list
km=modone((1:npts)-1,npts)';
k0=(1:npts)';
kp=modone((1:npts)+1,npts)';

for i=1:npts
    %get three points
    Xtemp=[Xcont(km(i)),Xcont(k0(i)),Xcont(kp(i))];
    Ytemp=[Ycont(km(i)),Ycont(k0(i)),Ycont(kp(i))];
    
    %find circle center (voronoi)
    try
        [Vx,Vy]=voronoi(Xtemp,Ytemp);
        
        %sig figs
        sig=5;
        Vx=reshape(round(Vx*10^sig)/10^sig,[],1);
        Vy=reshape(round(Vy*10^sig)/10^sig,[],1);
        
        compare_pts=zeros(6,6);
        %compare sets of points find intersection
        for j=1:6
            for k=1:6
                if and(and(Vx(j)==Vx(k),Vy(j)==Vy(k)),k~=j)
                    compare_pts(j,k)=1;
                end
            end
        end
        
        [ind1,ind2]=find(sum(compare_pts,1)==2,1,'first');
        VorX=Vx(ind1,ind2);
        VorY=Vy(ind1,ind2);
        
        %get circle radius
        R=sqrt((Xtemp(2)-VorX).^2+(Ytemp(2)-VorY).^2);
        
        %record curvature along contour
        kappa1(i)=1/R;
        
        %find normal axis vector
        norm1=[VorX-Xtemp(2),VorY-Ytemp(2)];
        norm1=norm1/norm(norm1)*maxL;
        
        %extend vector across centerline
        Xtend=[Xtemp(2)+norm1(1),Xtemp(2),Xtemp(2)-norm1(1)];
        Ytend=[Ytemp(2)+norm1(2),Ytemp(2),Ytemp(2)-norm1(2)];
        
        %find intersection of normal axis vector and centerline
        [Xint,Yint]=curveintersect(Xcent,Ycent,Xtend,Ytend);
        
        %record curvature
        if isempty(Xint)
            %no intersection
            kappa2(i)=NaN;
        else
            %find distances to intersection
            D=sqrt((Xtemp(2)-Xint).^2+(Ytemp(2)-Yint).^2);
            cond1=and(D>(1/2*1/max(abs(kappa1_other))),D<(1/2*maxL));
            cond2=min(D(cond1));
            
            if ~isempty(cond2)
                ind1=find(D==min(D(cond1)));
                
                %record the perpandicular curvature
                kappa2(i)=1/sqrt((Xtemp(2)-Xint(ind1))^2+(Ytemp(2)-Yint(ind1))^2);
                
                %{
            %test plot
            hold off
            plot(Xcont,Ycont,'k')
            hold on
            plot(Xcent,Ycent,'b')
            plot(Xtemp,Ytemp,'ro')
            plot([Xtemp(2),Xint(ind1)],[Ytemp(2),Yint(ind1)],'m')
            axis equal tight
            plot(Xint(ind1),Yint(ind1),'rd')
            drawnow
                %}
            else
                %no intersections are above the cutoff
                kappa2(i)=NaN;
            end
        end
    catch
        kappa2(i)=NaN;
    end
end

%***************************************
%******FIX NaNs & DICONTINUITIES********
%***************************************
try
    span=20;
    %pad input y data
    kappa2_pad=[kappa2(end-span+1:end);kappa2;kappa2(1:span)];
    
    %padded input x data
    xdata_pad=1:length(kappa2_pad);
    
    %find NaN values
    vec1=isnan(kappa2_pad);
    
    %{
    %define interpolation span
    span=21;
    degree=5;
    
    %find values that lie well outside local mean (2 STDs)
    kappa2_movavg=smooth(xdata_pad(~vec1),kappa2_pad(~vec1),span,'sgolay',degree);
    kappa2_dev=kappa2_pad(~vec1)-kappa2_movavg;
    %}
    
    std_cutoff=1/8;  %<-- use this value to adjust which points get rejected as being too noisy from the kappa2 data
    kappa2_pad_diff=diff([kappa2_pad;kappa2_pad(1)]);
    kappa2_pad_outliers=abs(kappa2_pad_diff)>(std_cutoff*std(kappa2_pad_diff));
    
    %construct final set of X Y input points
    kappa2_pad_select=kappa2_pad(and(~vec1,~kappa2_pad_outliers));
    xdata_pad_select=xdata_pad(and(~vec1,~kappa2_pad_outliers));
    
    %local interpolation of selected input data
    kappa2_pad_interp=interp1(xdata_pad_select,kappa2_pad_select,xdata_pad,'spline');
    
    %crop kappa2 vector back to proper length
    kappa2=kappa2_pad_interp(span+1:end-span);
catch
    disp('NaN correction encountered a problem.')
end

%{
%test plot for kappa smoothing
figure
plot(kappa2_pad,'kd')
hold on
plot(kappa2_pad_interp,'r')
axis tight
%}

%MODONE    Modulus after division starting from 1.
function r = modone(x,y)
if y==0
    r = NaN;
else
    n = floor(x/y);
    r = x - n*y;
    r(r==0)=y;
end






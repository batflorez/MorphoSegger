%Tristan Ursell
%March 2012
%
%This function takes a set points defined by X and Y, that are presumed to
%form a closed contour and returns a uniquely-determined mesh and branching
%neutral axis diagram, appended to the input data structure:
%
% frame.object
%
%Repeat points in (X,Y) will be removed.
%
%In particular, this script will add to this structure two sub-fields:
%
% frame.object.mesh
% frame.object.branches
%
%The 'mesh' sub-field has the sub-fields:
%
% X_pts = the x points of the current mesh polygon
% Y_pts = the y points of the current mesh polygon
% width = mean of the long side of the mesh block for the current point.
%     Assuming local azimuthal symmetry this should correspond to one of
%     the two principal curvatures, the other being the curvature along 
%     the contour.
% area = the area of the current mesh polygon
% perim = the perimeter length of the current mesh polygon
% centX = the x coordinate of the current mesh polygon centroid
% centY = the y coordinate of the current mesh polygon centroid
% neighbors = indices of the neighboring mesh polygons
%
%The 'branches' sub-field has the sub-fields:
%
% Xpos = the x coordinates of the current neutral axis
% Ypos = the y coordinates of the current neutral axis
% degree = the degree of connection of each points in the current neutral
% axis. A degree = 1 indicates a branch end-point that does not connect to
% another branch. A degree = 3 or higher indicates a branch end-point that
% meets that number of other branch end-points. A degree = 2 indicates the
% points is located within a branch.
% neighbors_left = the indices of all other branches that connect at one of
% the current branch.
% neighbors_right = the indices of all other branches that connect at the
% other end of the current branch.
%
%
%Example 1:
%
% N = 150;
% theta = linspace(0,2*pi,N);
%
% R = 4+cos(3*theta)+1/2*cos(theta)+1/4*cos(theta+pi/3);
% Y = 100 + 50*R.*sin(theta);
% X = 100 + 50*R.*cos(theta);
% O1 = interparc(N,X,Y,'linear');
% X = O1(1:end-1,1)';
% Y = O1(1:end-1,2)';
% [mesh,branches]=contour2mesh(X,Y,'plot');
%
%Example 2:
%
% Y = 100 + (10+1/4*sin(5*theta)).*sin(theta);
% X = 100 + (50+sin(5*theta)).*cos(theta);
% O1 = interparc(N,X,Y,'linear');
% X = O1(1:end-1,1)';
% Y = O1(1:end-1,2)';
% [mesh,branches]=contour2mesh(X,Y,'plot');


function [mesh,branches]=contour2mesh_cl(X,Y,varargin)

%check inputs
N=length(X);
if N~=length(Y)
    error('Input vectors must be the same length.')
end

if N<4
    error('Input vectors must have at least four points.')
end

%fix column vs. row
if size(X,1)>1
    X=X';
    Y=Y';
end

[~,I,~]=unique([X',Y'],'rows');
if length(I)~=N
    error('There can not be any repeat values in (X,Y).')
end

if size(X,1)>1
    X=X';
    Y=Y';
end

%{
%****** CELL WIDTH CALCULATIONS *****
%this analysis will fail on non-rod cells.

%set the minimum neighbor distance (i.e. do not connect neighbors below
%this distance)
min_d=pi*200/80;

if N>(2*min_d))
    %calculate the neighbor distance
    mat1=ones(1,N)'*(1:N);
    temp1=mat1-mat1';
    temp2=abs(sin(temp1*2*pi/(2*N)));
    neighbor_mat=temp2>sin(min_d*2*pi/(2*N));
    
    %calculate the distance matrix (row -> point of interest, column -> distance)
    Xpos_mat=ones(1,N)'*X;
    Ypos_mat=ones(1,N)'*Y;
    
    D=sqrt((Xpos_mat-Xpos_mat').^2+(Ypos_mat-Ypos_mat').^2); %symmetric
    
    %upsample and find minima not in the neighborhood region
    ribs=zeros(2,2,N);
    for i=1:N
        ytemp=-D(i,:);
        xtemp=1:length(ytemp);
        [xpeaks]=peakfind(xtemp,ytemp);
        
        if isempty(xpeaks)
            peak1=i;
        else
            %choose peak
            peaktemp=zeros(1,N);
            peaktemp(xpeaks)=1;
            peak1=find(peaktemp.*neighbor_mat(i,:));
            
            if isempty(peak1)
                peak1=i;
            end
        end
        
        %discard minima within neighborhood
        x1=X(i);
        x2=X(peak1);
        y1=Y(i);
        y2=Y(peak1);
        ribs(1,1,i)=x1;
        ribs(1,2,i)=y1;
        ribs(2,1,i)=x2;
        ribs(2,2,i)=y2;
    end
    
    %test figure
    figure;
    hold on
    plot(X,Y,'k-')
    plot(X(1),Y(1),'g.')
    plot(X(end),Y(end),'r.')
    for i=1:N
        plot(reshape(ribs(:,1,i),1,2),reshape(ribs(:,2,i),1,2),'r-')
    end
    axis equal tight
    xlabel('X')
    ylabel('Y')
end
%************************************
%}

%****** CELL MESHING AND BRANCHING *****
%get voronoi vertices
[Vx,Vy]=voronoi(X,Y);

%sig figs
sig=4;
Vx=round(Vx*10^sig)/10^sig;
Vy=round(Vy*10^sig)/10^sig;

%get voronoi cells
[V,C]=voronoin([X',Y']);
V_inside=find(inpolygon(V(:,1),V(:,2),X,Y))';

V=round(V*10^sig)/10^sig;

%determine which points lie inside the polygon
pts_inside=inpolygon(Vx,Vy,X,Y);

%find the segements that lie entirely inside the polygon
segs_inside=find(sum(pts_inside,1)==2);
n_in=length(segs_inside);

%find the segments that straddles the polygon
segs_strad=find(sum(pts_inside,1)==1);
n_strad=length(segs_strad);

%*****************************
%calculate radial segments
stradX=zeros(2,n_strad);
stradY=zeros(2,n_strad);

Xint=[X,X(1)]; %row vectors
Yint=[Y,Y(1)]; %row vectors
%first row = inside contour
%second row = at edge
for i=1:n_strad
    %get segment index
    col1=segs_strad(i);
    
    %get inside point
    stradX(1,i)=Vx(pts_inside(:,col1),col1);
    stradY(1,i)=Vy(pts_inside(:,col1),col1);
    
    %get outside point
%     [stradX(2,i),stradY(2,i)]=intersections(Vx(:,col1)',Vy(:,col1)',Xint,Yint); 
end

% comparing all the vectors in one function call greatly increases speed
% need to use NaN's to define discrete line segments (instead of a
% connected curve)
% output is ordered by the second curve, so need to re-index
NaNarray = nan(1,n_strad);
VxL = reshape([Vx(:,segs_strad);NaNarray],1,3*n_strad);
VyL = reshape([Vy(:,segs_strad);NaNarray],1,3*n_strad);
[X1,Y1,I1]=polyxpoly(VxL,VyL,Xint,Yint);
inds = ceil(I1(:,1)/3);
stradX(2,inds) = X1;
stradY(2,inds) = Y1;

%figure out which straddlers go with which contour points
%only two points can go with each contour point
for i=1:length(X)
    D=sqrt((X(i)-stradX(2,:)).^2+(Y(i)-stradY(2,:)).^2);
    [~,ord]=sort(D,'ascend');
    
    %get the indices of the strads
    ind1=ord(1:2);
    
    mesh(i).X_pts=[X(i),reshape(stradX(:,ind1),1,4)];
    mesh(i).Y_pts=[Y(i),reshape(stradY(:,ind1),1,4)];
    
    %calculate width / principal curvature (stradX/Y gives the sides
    %between contour and branch)
    xtemp=stradX(:,ind1);
    ytemp=stradY(:,ind1);
    
    dx1=xtemp(1,1)-xtemp(2,1);
    dx2=xtemp(1,2)-xtemp(2,2);
    dy1=ytemp(1,1)-ytemp(2,1);
    dy2=ytemp(1,2)-ytemp(2,2);
    
    mesh(i).width=1/2*(sqrt(dx1^2+dy1^2)+sqrt(dx2^2+dy2^2));
    
    %find current cell points that lie inside the contour
    rows=intersect(C{i},V_inside);
    
    Xtemp=V(rows,1)';
    Ytemp=V(rows,2)';
    
    %update cell point list
    mesh(i).X_pts=[mesh(i).X_pts,Xtemp];
    mesh(i).Y_pts=[mesh(i).Y_pts,Ytemp];
    
    %only keep the unique mesh points and order them with convexhulling
    temp=unique([mesh(i).X_pts',mesh(i).Y_pts'],'rows');
    K=convhull(temp(:,1),temp(:,2));
    
    mesh(i).X_pts=temp(K,1);
    mesh(i).Y_pts=temp(K,2);
    
    %calculate areas and centroids (using the exogeneous 'polygeom.m')
    [attribs,~,~]=polygeom(mesh(i).X_pts,mesh(i).Y_pts);
    mesh(i).area=attribs(1);
    mesh(i).perim=attribs(4);
    mesh(i).centX=attribs(2);
    mesh(i).centY=attribs(3);
end

%find adjacent mesh pieces
n_mesh=length(mesh);
adj_mat=false(n_mesh,n_mesh);

%find number of points in each mesh polygon
n = zeros(n_mesh,1);
for j=1:n_mesh
    n(j) = numel(mesh(j).X_pts(1:end-1));
end

%create index list
cumn = [0; cumsum(n)];
X2 = zeros(cumn(end),1);
Y2 = zeros(cumn(end),1);
ind = cell(n_mesh,1);
ind2 = zeros(cumn(end),1);
for j=1:n_mesh
    X2(cumn(j)+1:cumn(j+1)) = mesh(j).X_pts(1:end-1);
    Y2(cumn(j)+1:cumn(j+1)) = mesh(j).Y_pts(1:end-1);
    ind{j}=cumn(j)+1:cumn(j+1);
    ind2(cumn(j)+1:cumn(j+1))=j;
end

%do cross correlaton of 'inpolygon' points
for i=1:n_mesh
    %get points, arrange in matrix and remove endpoint
    X1=mesh(i).X_pts(1:end-1);
    Y1=mesh(i).Y_pts(1:end-1);
    
    %number of polygon sides (two sides are always on the outside)
    n_sides=length(X1)-2;
    
    %find neighbors
    [~,edge_list]=inpolygon(X2,Y2,X1,Y1);
    qs=0;
    [n_share,~] = histc(ind2(edge_list),1:n_mesh);
    adj_mat(i,setdiff(find(n_share>=2),i))=1;
%     for j=i:n_mesh
%         if i~=j
%             %if the points associated with the i-th and j-th share more
%             %than two points they are neighbors
%             if sum(edge_list(ind{j}))>=2
%                 qs=qs+1;
%                 adj_mat(i,j)=1;
%                 adj_mat(j,i)=1;
%             end
%             %if all the sides have been used up, exit the loop
%             if qs==n_sides
%                 break
%             end
%         end
%     end

    %record neighbors
    mesh(i).neighbors=find(adj_mat(i,:));
end

%************************************
%   PARSE NEUTRAL AXES
%************************************

%get internal segments
brch_X=Vx(:,segs_inside);
brch_Y=Vy(:,segs_inside);

%calculate number of distinct internal points and point connection matrix
points0=unique([reshape(brch_X,[],1),reshape(brch_Y,[],1)],'rows');
npts=size(points0,1);
point_mat=zeros(npts,npts);
for i=1:npts
    %get current point
    curr_pt=points0(i,:);
    
    %generate list of segments in voronoi space that share this point
    listx=brch_X==curr_pt(1);
    listy=brch_Y==curr_pt(2);
    list_both=listx.*listy;
    
    %indices of shared segments
    seg_list=find(sum(list_both,1)>0);
    
    %list of connected points in X Y space
    for j=1:length(seg_list)
        %ignore zero-length segments created by round-off error
        if ~and(brch_X(1,seg_list(j))==brch_X(2,seg_list(j)),brch_Y(1,seg_list(j))==brch_Y(2,seg_list(j)))
            %get distinct row
            row1=find(list_both(:,seg_list(j))==0);
            
            %record connected points
            xconn=brch_X(row1,seg_list(j));
            yconn=brch_Y(row1,seg_list(j));
            
            %find this point in the condensed list
            listx2=points0(:,1)==xconn;
            listy2=points0(:,2)==yconn;
            listboth2=find((listx2.*listy2)');
            
            %get list of points that are connected to current point
            point_mat(i,listboth2)=1;
        end
    end
end

%get point connection degrees
%degrees > 2 indicate a branch point
degrees=sum(point_mat,1);

if max(degrees(:))>3
    warning('Analysis found a branch point with more than three connections.')
end

%first start by handling end points of the network
ends=find(degrees==1)';
nodes=find(degrees>2);
node_num=length(nodes);

%create first grouping matrix (remove nodes from the regular diagram)
group_mat=point_mat+eye(npts);
group_mat(nodes',:)=0;
group_mat(:,nodes)=0;

%find matrix representation of nodes
node_mat=zeros(npts,npts);
node_mat(nodes',:)=1;
node_mat(:,nodes)=1;
node_mat=point_mat.*node_mat;

%find ends connected directly to nodes
lones=find(degrees==1); %these are the ends

q=0;
for i=1:length(lones)
    curr_vec=zeros(npts,1);
    curr_vec(lones(i))=1;
    
    end1=find(node_mat*curr_vec);
    
    if ~isempty(end1)
        q=q+1;
        branches(q).degree=[degrees(lones(i)),degrees(end1)];
        branches(q).Xpos=[points0(lones(i),1),points0(end1,1)];
        branches(q).Ypos=[points0(lones(i),2),points0(end1,2)];
    end
end

%find other groups
[groups0,~]=graph_analysis_morphometrics(group_mat);

%find connecting nodes for these groups
for i=1:length(groups0)
    q=q+1;
    
    curr_vec=zeros(npts,1);
    curr_vec(groups0(i).elements)=1;
    
    %find elements with nodes in this branch group
    elements=find((node_mat*curr_vec+curr_vec)>0);
    nodes1=find((node_mat*curr_vec)>0);
    ends1=intersect(ends,groups0(i).elements);
    
    %get points in this branch
    X0=points0(elements,1);
    Y0=points0(elements,2);
    
    %choose starting point
    if ~isempty(ends1)
        P=find(elements==ends1(1));
    else
        P=find(elements==nodes1(1));
    end
    
    %get ordering of points into a coherent branch
    [Xout,Yout,I]=points2contour(X0,Y0,P,'cw');
    
    branches(q).degree=degrees(I);
    branches(q).Xpos=Xout;
    branches(q).Ypos=Yout;
end

%find nodes connected directly to other nodes
if node_num>1
    for i=1:node_num
        %get connections to first node
        list1=[nodes(i),find(point_mat(nodes(i),:))];
        for j=1:node_num
            if j>i
                %get connections to second node
                list2=[nodes(j),find(point_mat(nodes(j),:))];
                
                %find overlapping nodes
                overlaps=intersect(list2,list1);
                
                if length(overlaps)>2
                    disp(['Something is up with node-node branch determination at node: ' num2str(i)])
                end
                
                %record overlaps (should be two)
                if ~isempty(overlaps)
                    q=q+1;
                    
                    %handle node-node connection
                    if length(overlaps)==2
                        branches(q).degree=degrees(overlaps);
                        branches(q).Xpos=points0(overlaps,1);
                        branches(q).Ypos=points0(overlaps,2);
                    else %handle node-point-node connection
                        p1=[nodes(i),overlaps(1),nodes(j)];
                        
                        branches(q).degree=degrees(p1);
                        branches(q).Xpos=points0(p1,1);
                        branches(q).Ypos=points0(p1,2);
                    end
                end
            end
        end
    end
end

%figure out which branches are neighbors
branchN=length(branches);
if branchN>1
    branchL=zeros(branchN,branchN);
    branchR=zeros(branchN,branchN);
    for i=1:branchN
        %compare LHS of branch to other branches
        xLi=branches(i).Xpos(1);
        yLi=branches(i).Ypos(1);
        
        %compare RHS of branch to other branches
        xRi=branches(i).Xpos(end);
        yRi=branches(i).Ypos(end);
        
        for j=1:branchN
            if i~=j
                %compare LHS of branch to other branches
                xj=[branches(j).Xpos(1),branches(j).Xpos(end)];
                yj=[branches(j).Ypos(1),branches(j).Ypos(end)];
                
                %branch touching conditions
                cond1=and(xLi==xj(1),yLi==yj(1));
                cond2=and(xLi==xj(2),yLi==yj(2));
                cond3=and(xRi==xj(1),yRi==yj(1));
                cond4=and(xRi==xj(2),yRi==yj(2));
                
                if or(cond1,cond2)
                    branchL(i,j)=1;
                end
                
                if or(cond3,cond4)
                    branchR(i,j)=1;
                end
            end
        end
        
        %get list of neighbors
        branches(i).neighbors_left=find(branchL(i,:));
        branches(i).neighbors_right=find(branchR(i,:));
    end
end

%calculate branch lengths
for i=1:length(branches)
    Xcurr=branches(i).Xpos;
    Ycurr=branches(i).Ypos;
    
    D=sqrt((Xcurr(1:end-1)-Xcurr(2:end)).^2+(Ycurr(1:end-1)-Ycurr(2:end)).^2);
    
    branches(i).length=sum(D);
end

%*******************************
%  OUTPUT PLOTS
if ~isempty(varargin)
    if strcmp('plot',varargin{1})
        h1=figure;
        %axes(handles.axes1);
        hold off
        plot(X,Y,'k-')
        hold on
        plot(X,Y,'k.')
        plot(stradX,stradY,'k-')
        plot(brch_X,brch_Y,'k-')
        xlim([min(X),max(X)])
        ylim([min(Y),max(Y)])
        axis equal
        box on
        
        cmap=jet(length(branches));
        for b=1:length(branches)
            Xbr=branches(b).Xpos;
            Ybr=branches(b).Ypos;
            plot(Xbr,Ybr,'Linewidth',2,'color',cmap(b,:));
            text(mean(Xbr),mean(Ybr),num2str(b))
        end
        
        %{
    n=1;
    for i=1:length(mesh(n).neighbors)
        m=mesh(n).neighbors(i);
        plot(mesh(m).X_pts,mesh(m).Y_pts,'r-','Linewidth',2)
    end
    plot(mesh(n).X_pts,mesh(n).Y_pts,'g-','Linewidth',2)
        %}
    end
end











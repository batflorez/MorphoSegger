%Tristan Ursell
%May 2014

function [frameout]=fit_contours(basename,frame,fr_num,Nim,handles,varargin)


global testV

global f_kint
global f_Nit
global f_r_int
global f_curve_std

global f_resize

global f_contmin
global f_acute
global f_pixelprox

seedq=get(handles.checkbox_seed,'Value');

if testV
    act_frame=str2num(get(handles.edit_test_frame,'String'));
else
    act_frame=fr_num;
end

set(handles.text_process,'String',['Proximity testing and contour fitting frame ' num2str(act_frame) ' of ' num2str(Nim)]);
drawnow
%%{
if and(or(~seedq,and(seedq,fr_num==1)),~testV)
    %***********************************************************
    %Calculate which cells are proximal to each other
    %***********************************************************
    %create proximity strel
    strel_prox=strel('disk',f_pixelprox);
    
    %load bwlabel image
    Im1=imread([basename '_Gparent.tif'],fr_num);
    
    %create distance transform image
    %bw_dist1=bwdist(Im1>0)<=f_pixelprox;
    
    %get bounding boxes for bw-regions
    [Sy,Sx]=size(Im1);
    prop1=regionprops(Im1,'BoundingBox');
    
    for q=1:frame.num_objs
        %get current box
        box1=cat(1,prop1(q).BoundingBox);
        
        %find bounds within dilation
        x1=round(box1(1)-2*f_pixelprox);
        x2=floor(box1(1)+box1(3)+2*f_pixelprox);
        y1=round(box1(2)-2*f_pixelprox);
        y2=floor(box1(2)+box1(4)+2*f_pixelprox);
        inds_y=y1:y2;
        inds_x=x1:x2;
        condy=and(inds_y>=1,inds_y<=Sy);
        condx=and(inds_x>=1,inds_x<=Sx);
        
        %select new box size
        inds_x = inds_x(condx);
        inds_y = inds_y(condy);
        
        %get sub-image
        subI=Im1(inds_y,inds_x);
        sub0=subI==q;
        sub1=imdilate(sub0,strel_prox)>0;
        
        %get subimage with nearby pixels labeled
        if isa(Im1,'uint16')
            sub2=uint16(~sub0.*sub1).*subI;
        else
            sub2=uint8(~sub0.*sub1).*subI;
        end
        
        %put proximal cell bw_id's into frame structure
        frame.object(q).proxID=unique(sub2(sub2~=0))';
        
        %check for stop
        global stopV
        if stopV
            set(handles.text_process,'String',['Proximity testing stopped.'])
            frameout=0;
            return
        end
    end
end
%}
%***********************************************************************
%Perform contour fitting using current pixel list as the initial guess
%***********************************************************************
%load and smooth image
%Im0=mat2gray(relnoise(imread(fname,i),filt1_sz,1,'disk'));

fname=[basename '.tif'];

%handle resizing
if f_resize~=1
    Im0=mat2gray(imresize(imread(fname,act_frame),f_resize,'bicubic'));
else
    Im0=mat2gray(imread(fname,act_frame));
end

if or(~seedq,and(seedq,fr_num==1))
    %create comparative image
    Iseg=imread([basename '_Gparent.tif'],fr_num)>0;
    bck1=Im0.*~Iseg;
    I2(:,:,1)=bck1;
    I2(:,:,2)=bck1+Im0.*Iseg;
    I2(:,:,3)=bck1;
    
    if max(Iseg(:))==0
        set(handles.text_process,'String','No objects found in this frame.')
        return
    end
else
    I2=Im0;
end

%plotting
%set(0,'CurrentFigure',gcf)
axes(handles.axes1);
hold off
axes(handles.axes1);
imagesc(I2,'Parent',handles.axes1)
colormap(gray)
hold on
axis equal tight

%convolve gaussian across entire image to get smoothed force fields
%get force field from magnitude image
Imag=mat2gray(imread([basename '_Gsegt.tif'],fr_num));
[Sy,Sx]=size(Imag);

%Image force smoothing filter

%{
%smooth then gradient
if f_r_int>0
    conv1_sz=ceil(3*f_r_int);
    if conv1_sz<=1
        conv1_sz=1;
    elseif mod(conv1_sz,2)==0
        conv1_sz=conv1_sz+1;
    end
    conv1=fspecial('gaussian',[conv1_sz,conv1_sz],f_r_int);
    norm1=conv2(ones(size(Imag)),conv1,'same');
    
    Imag_conv = conv2(Imag,conv1,'same')./norm1;
else
    Imag_conv=Imag;
end

[gx_conv,gy_conv]=gradient(Imag_conv);
%}

%%{
%gradient, then smooth
[gradX,gradY]=gradient(Imag);

%Image force smoothing filter
if f_r_int>0
    conv1_sz=ceil(3*f_r_int);
    if conv1_sz<=1
        conv1_sz=1;
    elseif mod(conv1_sz,2)==0
        conv1_sz=conv1_sz+1;
    end
    conv1=fspecial('gaussian',[conv1_sz,conv1_sz],f_r_int);
    norm1=conv2(ones(size(Imag)),conv1,'same');
    
    gx_conv = conv2(gradX,conv1,'same')./norm1;
    gy_conv = conv2(gradY,conv1,'same')./norm1;
else
    gx_conv = gradX;
    gy_conv = gradY;
end
%}

% take meshgrid indices for later use
[N,M] = size(gx_conv);
[x_mesh,y_mesh] = meshgrid(1:M,1:N);

%create list of objects that have cellIDs and possible repeat values
cell_vec=zeros(1,frame.num_objs);
for q=1:frame.num_objs
    %skip if the cell does not have an ID
    if ~isfield(frame.object(q),'cellID')
        continue
    elseif isempty(frame.object(q).cellID)
        continue
    end
    
    cell_vec(q)=frame.object(q).cellID;
end

%detect cells in frame that have the same cellID
cell_mat=zeros(length(cell_vec),length(cell_vec));
for i=1:length(cell_vec)
    for j=1:length(cell_vec)
        if and(and(cell_vec(i)==cell_vec(j),cell_vec(i)~=0),i~=j)
            cell_mat(i,j)=1;
        end
    end
end

%if cells have the same cellID throw warning
if sum(cell_mat(:))>1
    warning(['Objects with redundant cell IDs detected.'])
end

%examine each cell
for q=1:frame.num_objs
    %skip if the cell does not have an ID
    if cell_vec(q)==0
        continue
    end
    
    %check for stop
    global stopV
    if stopV        
        set(handles.text_process,'String',['Contour fitting stopped.'])
        frameout=0;
        return
    end
    
    clearvars X Y kappa peak1 peak2
    
    %initialize contour guess
    if and(seedq,fr_num>1)
        frame_m1=varargin{1};
        X=frame_m1.object(q).Xcont;
        Y=frame_m1.object(q).Ycont;
        
        %skip if contour doesn't exist
        if isempty(X)
            continue
        end
        
        X(end+1)=X(1);
        Y(end+1)=Y(1);
    else
        X=frame.object(q).Xperim;
        Y=frame.object(q).Yperim;
        
        %skip if contour doesn't exist
        if isempty(X)
            continue
        end
    end
    
    for k=1:f_Nit
        
        %calculate tangent lengths
        tvec=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
        contL=round(sum(tvec));
        
        %skip cell if contour gets too short
        if contL<(f_contmin*f_resize)
            set(handles.text_process,'String',['In frame ' num2str(act_frame) ', object ' num2str(q) ', the contour grew too short to process.'])
            break
        end
        
        %realign points
        try
            [X,Y]=interpcurv2n(X,Y,contL);
        catch
            O1 = interparc(contL,X,Y,'linear');
            X=O1(:,1);
            Y=O1(:,2);
        end
        
        %take care of points on boundary
        X(X<1)=1;
        Y(Y<1)=1;
        X(X>Sx)=Sx;
        Y(Y>Sy)=Sy;
        
        %get box size for current cell
        inds_x = round(min(X)-f_pixelprox):round(max(X)+f_pixelprox);
        inds_y = round(min(Y)-f_pixelprox):round(max(Y)+f_pixelprox);
        
        %correct if box is on border
        condy=and(inds_y>=1,inds_y<=Sy);
        condx=and(inds_x>=1,inds_x<=Sx);
        
        %select new box size
        inds_x = inds_x(condx);
        inds_y = inds_y(condy);
        
        %image forces calculations
        %get sub mesh from gradient for current cell
        subX=x_mesh(inds_y,inds_x);
        subY=y_mesh(inds_y,inds_x);
        
        try
            %calculate image forces from the interpolated grid
            Fx_int = interp2a(subX,subY,gx_conv(inds_y,inds_x),X,Y);
            Fy_int = interp2a(subX,subY,gy_conv(inds_y,inds_x),X,Y);
            
            %adjust point positions
            X=X+f_kint*Fx_int;
            Y=Y+f_kint*Fy_int;
        catch
            display(['warning:  force field interpolation error, frame ' num2str(act_frame) ', object ' num2str(q)]);
        end
        
        %stochastic noise-based annealing of contour
        %X=X+f_kint*Fx_int+normrnd(0,2*exp(-5*(k-1)/f_Nit)*f_kint,size(Fx_int));
        %Y=Y+f_kint*Fy_int+normrnd(0,2*exp(-5*(k-1)/f_Nit)*f_kint,size(Fy_int));
    end
    
    %skip cell if contour gets too short
    if contL<(f_contmin*f_resize)
        continue
    end
    
    %display and update title
    title(['Frame ' num2str(act_frame) ' of ' num2str(Nim) ', cell ID: ' num2str(frame.object(q).cellID)])
    
    %calculate tangent lengths
    tvec=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
    
    %realign points (currently wrapped contour)
    try
        [X,Y]=interpcurv2n(X/f_resize,Y/f_resize,round(sum(tvec)/f_resize));
    catch
        O1 = interparc(round(sum(tvec)/f_resize),X/f_resize,Y/f_resize,'linear');
        X=O1(:,1);
        Y=O1(:,2);
    end
    
    %calculate tangent vector lengths for wrapped contour
    tmag=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
    
    %impose minimum contour length
    if length(X)>=f_contmin
        %enter loop to get rid of kinks in contour
        q1=0;
        while length(X)-q1>(2*f_contmin)
            q1=q1+1;
            
            %set precision
            ndigs=8;
            X=round(X*10^ndigs)/10^ndigs;
            Y=round(Y*10^ndigs)/10^ndigs;
            
            %calculate raw contour curvature (unwrapped)
            kappa0=contour2curvature(X,Y);
            
            %list of points to be included (those where tangents are not at
            %right angles) (s1 = 1 excludes acute angles, s1 > 1
            %excludes angles greater than 90')
            s1=f_acute; %<------------ maximum allowable curvature
            inc_list=1./abs(kappa0)/mean(tmag)>s1/sqrt(2);
            
            if sum(inc_list)==length(inc_list)
                break
            else
                %exclude points based on curvature
                pt1=find(inc_list,1,'first');
                X=[X(inc_list);X(pt1)];
                Y=[Y(inc_list);Y(pt1)];
                
                %calculate tangent vector lengths for wrapped contour
                tmag=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
                
                %realign points
                if length(X)<f_contmin
                    break
                else
                    try
                        [X,Y]=interpcurv2n(X,Y,round(sum(tmag)));
                    catch
                        O1 = interparc(round(sum(tmag)),X,Y,'linear');
                        X=O1(:,1);
                        Y=O1(:,2);
                    end
                end
            end
        end
        
        %edge case that original contour is a line that collapses to a
        %point when acute angles removed
        if isempty(X)
            continue
        end
        
        %unwrap contour
        X=X(1:end-1);
        Y=Y(1:end-1);
        
        %remove any duplicate points
        [~,I,~]=unique([X,Y],'rows','stable');
        if length(I)~=length(X)
            warning('Detected repeat values in (X,Y).')
            X=X(I);
            Y=Y(I);
        end
        
        %wrap contour
        X=[X;X(1)];
        Y=[Y;Y(1)];
        
        %calculate new raw contour curvature
        kappa=contour2curvature(X,Y);
        
        %calculate contour length using wrapped contour
        tmag=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
        contL=sum(tmag); %total arc-length
        frame.object(q).arc_length=tmag; %differential arc-length
        frame.object(q).area=polyarea(X,Y); %in pixels^2
        frame.object(q).circularity=contL^2/frame.object(q).area/(4*pi); %circularity

        %unwrap contour
        X=X(1:end-1);
        Y=Y(1:end-1);
        
        %pole finding and curvature smoothing <------------
        if length(X)>(5*f_contmin)
            kappaX=1:3*length(kappa);
            kappaY=[kappa;kappa;kappa]';
            pole_sz=round(4*f_contmin);
            if mod(pole_sz,2)==0
                pole_sz=pole_sz+1;
            end
            [~,yout,peakspos]=peakfind(kappaX,kappaY,1,pole_sz,f_contmin);
            
            %find two highest peak indices
            yout=yout(length(kappa)+1:2*length(kappa));
            peakspos=peakspos(and(peakspos>=(length(kappa)+1),peakspos<=(2*length(kappa))))-length(kappa);
            peakheights=yout(peakspos);
            temp1=sort(peakheights,'descend');
            
            if length(peakspos)>1
                %check how many peaks meet the first and second height
                peak_inds1=find(yout==temp1(1));
                peak_inds2=find(yout==temp1(2));
                
                if length(peak_inds1)>1
                    peak1=peak_inds1(1);
                    peak2=peak_inds1(2);
                else
                    peak1=peak_inds1(1);
                    peak2=peak_inds2(1);
                end
         
                %get contour positions of peaks
                x_peak_1=X(peak1);
                y_peak_1=Y(peak1);
                x_peak_2=X(peak2);
                y_peak_2=Y(peak2);             
                
                dpX=abs(x_peak_1-x_peak_2);
                dpY=abs(y_peak_1-y_peak_2);
                
                p1=sqrt(x_peak_1^2+y_peak_1^2);
                p2=sqrt(x_peak_2^2+y_peak_2^2);
                
                %pick peak based on relative position
                if abs(p1-p2)<(2*f_contmin)
                    if dpX<20 %if the X positions are too close, choose
                        %based on farthest Y value
                        if y_peak_1>y_peak_2
                            ind_rot=peak1;
                        else
                            ind_rot=peak2;
                        end
                    else %if the X positions are far enough, choose one
                        %with X position farthest from origin
                        if x_peak_1>x_peak_2
                            ind_rot=peak1;
                        else
                            ind_rot=peak2;
                        end
                    end
                else
                    if p1>p2 %pick the one farthest from the image origin
                        ind_rot=peak1;
                    else
                        ind_rot=peak2;
                    end
                end
            elseif length(peakspos)==1 %if there is only one peak, use it
                ind_rot=find(yout==temp1(1),1,'first');
                peak1=ind_rot;
            else
                ind_rot=0;
            end
            
            %curvature smoothing
            curve_sz=ceil(2*3*f_curve_std);
            if mod(curve_sz,2)==0
                curve_sz=curve_sz+1;
            end
            [~,yout2,~]=peakfind(kappaX,kappaY,1,curve_sz,f_curve_std);  %<------ std of gaussian smoothing filter, in pixels
            
            %rotate contour
            Xrot=circshift(X,[-(ind_rot),0]);
            Yrot=circshift(Y,[-(ind_rot),0]);
            kappa_rot=circshift(yout2(length(kappa)+1:2*length(kappa))',[-(ind_rot-1),0]);
            kappa_raw=circshift(kappa,[-(ind_rot-1),0]);
            
            %assign poles indices
            if length(peakspos)>1
                frame.object(q).pole1=mod(peak1-ind_rot,length(Xrot))+1;
                frame.object(q).pole2=mod(peak2-ind_rot,length(Xrot))+1;
            elseif length(peakspos)==1
                frame.object(q).pole1=mod(peak1-ind_rot,length(Xrot))+1;
                frame.object(q).pole2=round(mod(peak1-ind_rot,length(Xrot))+1+contL/2);
            else
                frame.object(q).pole1=1;
                frame.object(q).pole2=round(contL/2);
            end
        else
            %if the contour is too short, don't align or smooth
            Xrot=X;
            Yrot=Y;
            kappa_rot=kappa;
            kappa_raw=kappa;
            frame.object(q).pole1=1;
            frame.object(q).pole2=round(contL/2);
        end
        
        %update data structure
        frame.object(q).kappa_raw=kappa_raw;
        frame.object(q).kappa_smooth=kappa_rot;
        frame.object(q).Xcont=Xrot;
        frame.object(q).Ycont=Yrot;
        frame.object(q).num_pts=length(Xrot);
        
        %calculate centroid
        frame.object(q).Xcent_cont=mean(Xrot);
        frame.object(q).Ycent_cont=mean(Yrot);
        
        %calculate angle via Deming Regression
        n=length(Xrot);
        xcent=Xrot-mean(Xrot);
        ycent=Yrot-mean(Yrot);
        s_xx=1/(n-1)*sum((xcent-mean(xcent)).^2);
        s_yy=1/(n-1)*sum((ycent-mean(ycent)).*(xcent-mean(xcent)));
        s_xy=1/(n-1)*sum((xcent-mean(xcent)).^2);
        delta=1;
        
        m=(s_yy-delta*s_xx+sqrt((s_yy-delta*s_xx)^2+4*delta*s_xy^2))/(2*s_xy);
        frame.object(q).theta_cont=atan2(m,1);
        
        %graphics
        plot(Xrot*f_resize,Yrot*f_resize,'c-','LineWidth',1,'Parent',handles.axes1)
        plot(Xrot(end)*f_resize,Yrot(end)*f_resize,'r.','Parent',handles.axes1)
        plot(Xrot(1)*f_resize,Yrot(1)*f_resize,'g.','Parent',handles.axes1)
        if frame.object(q).pole1==1
            plot(Xrot(frame.object(q).pole2)*f_resize,Yrot(frame.object(q).pole2)*f_resize,'rx','Parent',handles.axes1)
        else
            plot(Xrot(frame.object(q).pole1)*f_resize,Yrot(frame.object(q).pole1)*f_resize,'rx','Parent',handles.axes1)
        end
        text(mean(Xrot)*f_resize,mean(Yrot)*f_resize,num2str(frame.object(q).cellID),'color',[1,1,1],'Parent',handles.axes1)
        drawnow
    else
        frame.object(q).Xcont=X;
        frame.object(q).Ycont=Y;
        frame.object(q).num_pts=length(X);
        
        %calculate centroid
        frame.object(q).Xcent_cont=mean(X);
        frame.object(q).Ycent_cont=mean(Y);
    end
end

frameout=frame;








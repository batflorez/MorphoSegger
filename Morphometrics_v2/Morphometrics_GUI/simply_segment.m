function [segdata,Imag_out,Iseg_out,Im1,Icolor]=simply_segment(Im0,mask1,fr_num,rect1,handles)
%Tristan Ursell
%April 2013
%
%This m-file peforms gradient segmentation of images presumed to be phase
%dark -- to use phase light or internal fluorescence, simply invert the
%image. The more cells in each image the longer it will take.  The larger
%the image, regardless of object number, the longer it will take.
%
% Im0 is an input image matrix of any size or class.
%
% Image parameters:
% 'resize' = scaling factor for image (0 .. inf), default is 1.
% 'histlim' = percentage of histogram to saturate on either side (0 .. 0.5)
% default is 0.001.
% 'filtsz' = size of relative noise filter to apply, must be odd, defualt
% is 3.
%
% Segmentation parameters:
% 'hmin' = h-min transform parameter (0..1) sets depth threshold for
% watershed segementation, default is 0.1.
% 'hmin_split' = secondary h-min transform parameter (0..hmin), default is
% off.
% 'rejfalpos' = 1 rejects false positives (according to internal criteria),
% default is off.
%
% Object parameters:
% 'areamin' = minimum object area in pixels, default is 100.
% 'areamax' = maximum object area in pixels, default is image size / 100.
%
% Output options:
% 'plot_out' will display the segmented image.
% 'im_out' will save the logical segmented image.
%
%The structure array 'segdata' contains:
%
%segdata.num_cells = number of identified objects.
%
%segdata.Im_out = segmented image.
%
%segdata.cell(i).Xcent = x coordinate of object 'i' centroid.
%segdata.cell(i).Ycent = y coordinate of object 'i' centroid.
%
%segdata.cell(i).area = area of object 'i' in pixels.
%
%segdata.cell(i).theta = object orientation in radians with respect to
%image coordiantes.
%
%

%**************************************************************************
global testV

%type of segmentation (1 = gradient, 2 = laplacian, 3 = adaptive threshold)
global v_method

%type of image (1 = phase dark, 2 = interior fluor, 3 = peripheral fluor)
global v_imtype

%put contrast in the correct direction for processing
if v_imtype==2
    Im0=mat2gray(-single(Im0));
else
    Im0=mat2gray(single(Im0));
end

%------ SEGEMENTATION------
%image scaling / resizing factor
global f_resize

%background removal filter size (pixels)
global f_back
if f_back>0
    if and(f_back<=3,f_back>0)
        f_back=3;
    elseif and(mod(f_back,2)==0,f_back>0)
        f_back=f_back+1;
    end
    set(handles.edit_back,'String',num2str(f_back))
end

%size of the relnoise disk filter
global f_relsz

%size of the relnoise disk filter
global f_int_rej

%minimum relative variance required for any segmentation
global f_seg_min

%histogram saturation percentage
global f_histlim

%minimum undilated cell area
%f_areamin=200;
global f_areamin

%maximum undilated cell area
%f_areamax=50*f_areamin;
global f_areamax

%'hmin' transform thresholds
global f_hmin
global f_hmin_split

%gradient smoothing filter
%f_gstd=1;
global f_gstd

%perimeter dilation for including gradient points for fitting (pixels)
if f_gstd<2
    f_dilt=5;
else
    f_dilt=ceil(3*f_gstd);
end

%is this a seeded contour finding?
seedq=get(handles.checkbox_seed,'Value');

%**************************************************************************
%determine if the image has sufficient intensity difference between
%median intensity and min and max intensity to warrant segmentation

%{  
%old way
%calculate median intensity of image
med_int=median(Im0(:));
lim2=stretchlim(Im0,[0.001 0.999]);
%lim2=stretchlim(Im0,[0.005 0.995]);
max_int=mean(Im0(Im0(:)>=lim2(2)));
min_int=mean(Im0(Im0(:)<=lim2(1)));

ratio1=max_int/med_int;
ratio2=med_int/min_int;
%}

%calculate median intensity of image
lim2=stretchlim(Im0,[f_seg_min 1-f_seg_min]);

if v_imtype==1
    reg_test=imfill(Im0<=lim2(1),'holes');
else
    reg_test=imfill(Im0>=lim2(2),'holes');
end

reg_test_prop=regionprops(reg_test,'Area');

%check to make sure there is enough variation in the image to warrant segmetation
if max([reg_test_prop.Area])<f_areamin %<---- sets the discrimination for segmentation
    Imag_out=uint16(zeros(size(Im0)));
    Iseg_out=uint16(zeros(size(Im0)));
    Im1=Im0;
    
    segdata.num_objs=0;
    set(handles.text_process,'String','Image did not meet the minimum requirements for segmentation.')
    if testV
        pause(2)
    end
    
    Icolor=Im0;

    %reverse contrast when appropriate
    if v_imtype==2
        Icolor=1-Icolor;
    end
    
    %record output images
    Iseg_out=uint16(zeros(size(Im0)));
    Imag_out=uint16(zeros(size(Im0)));
    
    return
end
%**************************************************************************

%apply background removal (high pass gaussian)
if f_back>1
    set(handles.text_process,'String','Removing background ...')
    drawnow
    
    filt1=fspecial('gaussian',[round(3.5*f_back), round(3.5*f_back)],f_back); %<--- adjust relative size of guassian blur filter
    
    whitenorm1=conv2(ones(size(Im0)),filt1,'same');
    Im_bck=conv2(Im0,filt1,'same')./whitenorm1;
    
    Im0=Im0-Im_bck;
    
    if any(Im0(:)<0)
        disp('Background removal created some negative image values.')
    end
end

%apply gaussian smoothing
if f_gstd>0.5
    set(handles.text_process,'String','Smoothing image ...')
    drawnow
    
    gstd_sz=round(3.5*f_gstd);
    if mod(gstd_sz,2)==0
        gstd_sz=gstd_sz+1;
    end
    
    filt2=fspecial('gaussian',[gstd_sz,gstd_sz],f_gstd); %<--- adjust relative size of guassian blur filter
    
    whitenorm2=conv2(ones(size(Im0)),filt2,'same');
    Im0=conv2(Im0,filt2,'same')./whitenorm2;
end

%apply relnoise filter
[I1,Ivar,~]=relnoise(Im0,f_relsz,1,'disk');  %<--- consider using variance instead of grad magnitude?

%***************************************
%help the user decide on a primary threshold using the distribtuion of local image variances
if and(get(handles.checkbox_help,'Value'),v_method==1)
    %find distribution of local variances
    %the first peak of Nvar should approx. equal the primary threshold
    
    %find spatially varying image STD
    Istd=sqrt(Ivar);
    
    [Nvar,Xvar]=hist(Istd(:),round(min(size(I1))/2));
    
    h1=figure;
    semilogy(Xvar(Nvar>0),1+Nvar(Nvar>0),'k.','MarkerSize',12)
    axis tight
    xlabel('Local Standard Deviation (zero is saturated pixels)')
    ylabel('Number of Pixels')
    title('Estimate primary threshold from far Y=0 intercept (close to continue).')
    drawnow
    uiwait(h1)
    
    prompt={'Enter estimate of primary threshold:'};
    name='Primary Threshold Estimation';
    numlines=1;
    defaultanswer={num2str(round(2/3*100*Xvar(end-2))/100)};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answer)
        f_hmin=str2double(answer{1});
        set(handles.edit_primary,'String',num2str(f_hmin))
    else
        segdata=[];
        
        Icolor=Im0;
        
        %reverse contrast when appropriate
        if v_imtype==2
            Icolor=1-Icolor;
        end
        
        %record output images
        Iseg_out=uint16(zeros(size(Im0)));
        Imag_out=uint16(zeros(size(Im0)));
        
        return
    end
    
elseif and(get(handles.checkbox_help,'Value'),v_method==3)
    %extract sub-region
    temp2=I1(rect1(2):rect1(2)+rect1(4),rect1(1):rect1(1)+rect1(3));
    
    %calculate background means and stds
    bck_mean1=mean(temp2(:));
    bck_std1=std(double(temp2(:)));
    
    Itest1=(I1-bck_mean1)/bck_std1;
    
    [Nvar,Xvar]=hist(-Itest1(:),round(min(size(I1))/3));
    Xvar=Xvar(Nvar>0);
    Nvar=Nvar(Nvar>0);
    
    %find peaks
    peaks1=peakfind(Xvar,Nvar);
    
    %estimate threshold
    region1=and(Xvar>0,Xvar<=peaks1(end));
    sort1=Nvar(region1);
    mins1=sort(sort1,'ascend');
    posits1=region1.*any([mins1(1)==Nvar;mins1(2)==Nvar;mins1(3)==Nvar]);
    thres_est=mean(Xvar(posits1==1));
    
    h1=figure;
    semilogy(Xvar,1+Nvar,'k.','MarkerSize',12)
    hold on
    semilogy([0,0],[1,max(1+Nvar)],'r-')
    semilogy(Xvar(posits1==1),1+Nvar(posits1==1),'bx')
    semilogy([thres_est,thres_est],[1,max(1+Nvar)],'b-')
    axis tight
    xlabel('Local Standard Deviation (zero is saturated pixels)')
    ylabel('Number of Pixels')
    title('Estimate primary threshold from trough between the zero and positive peaks. (close to continue).')
    drawnow
    uiwait(h1)
    
    prompt={'Enter estimate of primary threshold:'};
    name='Primary Threshold Estimation';
    numlines=1;
    defaultanswer={num2str(round(thres_est*100)/100)};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answer)
        f_hmin=str2double(answer{1});
        set(handles.edit_primary,'String',num2str(f_hmin))
    else
        segdata=[];
        
        Icolor=Im0;
        
        %reverse contrast when appropriate
        if v_imtype==2
            Icolor=1-Icolor;
        end
        
        %record output images
        Iseg_out=uint16(zeros(size(Im0)));
        Imag_out=uint16(zeros(size(Im0)));
        
        return
    end
end

%hmin help only applies to gradient segmentation
if and(get(handles.checkbox_help,'Value'),v_method==2)
    set(handles.text_process,'String','warning: Threshold help does not apply to this method of segmentation.')
    set(handles.checkbox_help,'Value',0)
end
%***************************************

%resize image after initial manipulations so as not to introduce artifacts
if f_resize~=1
    set(handles.text_process,'String','Resizing image ...')
    drawnow
    
    I1=imresize(I1,f_resize,'bicubic');
    mask1=imresize(single(mask1),f_resize,'nearest');
end

%adjust object size contraints appropriately
areamin_use=round(f_resize^2*f_areamin);
areamax_use=round(f_resize^2*f_areamax);

%create pre-output magnitude image
if or(v_imtype==1,v_imtype==2)
    %calculate gradient
    [gx,gy]=gradient(I1);
    
    magG=sqrt(gx.^2+gy.^2);
else
    magG=I1;
end
%magG=relnoise(magG,f_relsz,1,'disk'); %possibly turn off

%apply contrasting and change class
if and(f_histlim>0,and(v_method~=2,v_method~=4))
    lim1=stretchlim(I1,[f_histlim, 1-f_histlim]);
    I1=imadjust(I1,[lim1(1),lim1(2)],[0 1]);
end

%use as output pre-process image
Im1=I1;

if or(~seedq,and(seedq,fr_num==1))
    %***************************************
    %        gradient pre-filtering
    %***************************************
    if v_method==1
        %create pre-output magnitude image
        if or(v_imtype==1,v_imtype==2)
            %calculate gradient
            [gx,gy]=gradient(I1);
            mag_temp=sqrt(gx.^2+gy.^2);
        else
            mag_temp=I1;
        end
        %magG=relnoise(magG,f_relsz,1,'disk'); %possibly turn off
        
        %apply hmin transform for watershedding to smoothed input image
        G_min=imhmin(mag_temp,f_hmin,8);
        
        %watershed hmin input image
        %G_shed=watershed(G_min,4);
        G_shed=watershed(G_min,8);
        
        %apply area constraints to watershed
        G_seg=bwareaopen(G_shed,areamin_use,4);
        %G_seg=bwareaopen(G_shed,areamin_use,4)&~bwareaopen(G_shed,areamax_use,4);
        
        %apply mask and round off and get rid of spurs
        G_fill=bwmorph(G_seg,'spur','Inf').*mask1;
        %G_fill=imfill(bwmorph(G_seg,'spur','Inf'),'holes');
        
        %figure;imagesc((G_fill>0).*G_min)
    end
    %---------------------------------------
    %---------------------------------------
    
    
    %***************************************
    %       laplacian pre-filtering
    %***************************************
    if v_method==2
        %take Laplacian of input image, normalize, and smooth
        L0=zerocross(del2(I1),1);
        
        %remove bits from inside of cells
        I2=bwareaopen(L0,areamin_use,8);
        
        %remove spur and h-bridged pixels
        I3=bwmorph(~bwmorph(I2,'spur','Inf'),'spur','Inf');
        G_shed=bwmorph(~bwmorph(~I3,'majority'),'spur','Inf');
        
        %remove small areas and big areas, clear border
        G_seg=bwareaopen(G_shed,areamin_use,4);
        %I5=bwareaopen(I4,areamin_use,4)&~bwareaopen(I4,areamax_use,4);
        
        %fill holes below a certain size and break cells free from background
        da=areamin_use;
        G_fill=bwareaopen(imopen(~bwareaopen(~G_seg,da,4),strel('disk',2)),areamin_use,4).*mask1;
    end
    %---------------------------------------
    %---------------------------------------
    
    
    %***************************************
    %   adaptive threshold pre-filtering
    %***************************************
    if v_method==3
        if get(handles.checkbox_simplethres,'Value')
            G_shed=I1<=f_hmin;
        else
            %extract sub-region
            temp2=I1(rect1(2):rect1(2)+rect1(4),rect1(1):rect1(1)+rect1(3));
            
            %calculate background means and stds
            bck_mean1=mean(temp2(:));
            bck_std1=std(double(temp2(:)));
            
            %calculate STD thresheld image
            %G_fill=(bck_mean1-I1)>(f_hmin*bck_std1);
            G_shed=(bck_mean1-I1)>(f_hmin*bck_std1);
        end
        
        if v_imtype==3
            G_shed=imfill(G_shed,'holes',4);
        end
        
        %apply area constraints to thresheld image
        G_seg=bwareaopen(G_shed,areamin_use,4);
        %G_seg=bwareaopen(I_thres,areamin_use,4)&~bwareaopen(I_thres,areamax_use,4);
        
        %apply mask and fills holes
        G_fill=imfill(G_seg,4,'holes').*mask1;
    end
    %---------------------------------------
    %---------------------------------------
    
    %***************************************
    %   canny edge pre-filtering
    %***************************************
    if v_method==4
        global f_canny_th1
        global f_canny_th2
        global f_canny_sig
        
        %find edges
        I2=edge(I1,'canny',[f_canny_th1,f_canny_th2],f_canny_sig);
        
        %remove spurs
        I3=bwmorph(I2,'spur','Inf');
        
        %apply min area constraint and fill holes
        G_shed=imfill(I3,'holes');
        G_seg=bwareaopen(G_shed,areamin_use,4);
        G_fill=G_seg;
    end
    %---------------------------------------
    %---------------------------------------
    
    %get rid of very large false positive regions (e.g. for topological
    %reasons, sometimes the background is segmented as its own region
    Ltemp=bwlabel(bwareaopen(G_fill,4*areamax_use,8));
    if max(Ltemp(:))==1
        G_fill=G_fill.*(~Ltemp);
    end
    
    %*****************************************
    %     SECONDARY SEGMENTATION
    %*****************************************
    %>>>> NOTES:  made this based on min pixel approach between sides %<<<<<<
    if and(f_hmin_split>=1,get(handles.checkbox_prox,'Value'))
        
        %get distance transform from region borders
        Gdist=-bwdist(~G_fill);
        
        %segment current image
        %G_fill2=G_fill.*(watershed(imhmin(Gdist,f_hmin_split),4)>0);
        Gdist(Gdist<=-f_hmin_split)=-f_hmin_split;
        G_fill2=G_fill.*(watershed(Gdist,4)>0);
        
        %get mask of regions that are too small to be sub-divided
        small_mask=bwlabel(G_fill - bwareaopen(G_fill,2*areamin_use,8));
     
        %put back regions that were too small to be resegmented
        for j=1:max(small_mask(:))
            cond1=G_fill2.*(small_mask==j);
            if sum(cond1(:))==0
                G_fill2=G_fill2 + (small_mask==j);
            end
        end
    else
        G_fill2=G_fill;
    end
    %*****************************************
    %*****************************************
    
    %{
    %clear encompassing border regions greater than a certain size
    G_fill2=imclearborder(G_fill2,4);
    G_bord=bwlabel(G_fill2);
    temp1=regionprops(logical(G_bord),'Area');
    bord_area=find(([temp1.Area]/(size(G_bord,1)*size(G_bord,2)))>0.25);
    
    %set large border regions to zero
    for j=1:length(bord_area)
        G_fill2(G_bord==bord_area(j))=0;
    end
    %}
    
    %clear border and larger interior regions
    if get(handles.checkbox_excludeedge,'Value')
        G_fill2=imclearborder(G_fill2.*(~bwareaopen(G_fill2,areamax_use,4)),4);
    else
        G_fill2=G_fill2.*(~bwareaopen(G_fill2,areamax_use,4));
    end
    %G_fill2=imfill(imclearborder(G_fill2,4),8,'holes');
    %figure; imagesc(G_fill+G_fill2)
    
    if max(G_fill2(:))>0
        %remove false positivess
        %int_rej sets a threshold for false positive detection
        %int_rej < 0 will encourage false positives.
        %0 < int_rej < 1 will weakly reject false positives
        %int_rej = 1 should be a decent choice to reject actual false positives
        %int_rej > 1 will tend to create false negatives
        
        if get(handles.checkbox_falsepos,'Value')
            if v_method==3
                if ~get(handles.checkbox_simplethres,'Value')
                    Lout=falsepos(G_fill2,I1,f_int_rej,bck_std1);
                else
                    Lout=falsepos(G_fill2,I1,f_int_rej,f_hmin);
                end
            elseif or(v_method==1,v_method==2)
                Lout=falsepos(G_fill2,I1,f_int_rej,f_hmin);
            else
                Lout=G_fill2;
            end
        else
            Lout=G_fill2;
        end
        
        if max(Lout(:))==0
            set(handles.text_process,'String',['False positive rejection is removing all elements, consider deselecting this feature.'])
            if testV
                pause(2)
            end
        end
        
        %expand the borders on segmented objects when peripheral fluor
        if and(v_method==2,v_imtype==3)
            Lout=bwmorph(Lout,'thicken',3);
        end
    else
        Lout=zeros(size(I1));
    end
    
    %create final output gradient magnitude image
    magG_out=mat2gray(imfill(imdilate(Lout,strel('disk',7)),'holes').*magG); %<---- adjust dilation around magG image
    
    %final removal of out of size contraint regions
    L_final=bwlabel(bwareaopen(Lout,areamin_use,4)&~bwareaopen(Lout,areamax_use,4),4);
    
    %construct color-coded output image
    %I1 = input image after various manipulations, pre-segmentation
    %G_shed = the raw segmented image
    %G_seg = the image after minimum size removal
    %G_fill2 = the image after maximum size removal and border clearing, but
    %before false-positive removal
    %L_final = final segmented image (that moves onto the next stage)
    true1=L_final>0;
    
    %creat output plot
    mask2=double(mask1);
    mask2(mask1==0)=0.5;
    mask2(1,1)=0;
    
    if get(handles.checkbox_colorcode,'Value')
        %false1=(G_fill2>0).*(~true1);
        false1=(G_fill2>0).*(~Lout);
        
        switch v_method
            case 1
                size1=~(G_shed<=1).*(~true1).*(~false1);
            case 2
                size1=(G_shed==0).*(~true1).*(~false1);
            case 3
                size1=(G_shed>0).*(~true1).*(~false1);
            case 4
                size1=~(G_shed<=1).*(~true1).*(~false1);
        end
        
        bck1=mat2gray(I1).*(~size1).*(~false1).*(~true1);
        Icolor(:,:,1)=(bck1+size1+false1).*mask2;
        Icolor(:,:,2)=(bck1+false1+true1).*mask2;
        Icolor(:,:,3)=(bck1).*mask2;
    else
        bck1=mat2gray(I1).*(~true1);
        Icolor(:,:,1)=(bck1).*mask2;
        Icolor(:,:,2)=(bck1+true1).*mask2;
        Icolor(:,:,3)=(bck1).*mask2;
    end
else
    Icolor=I1;
    L_final=zeros(size(Icolor));
    magG_out=mat2gray(magG).*mask1;
end

%reverse contrast when appropriate
if v_imtype==2
    Icolor=1-Icolor;
end

%record output images
Iseg_out=uint16(L_final);
Imag_out=im2uint16(magG_out);

%***************** HERE IS WHERE IT USES THE MASK ********************
%From Here I can use this to fit contours

if max(L_final(:))>0
    %find number of regions and region properties
    Ctemp=regionprops(L_final,'Centroid','Area','Orientation');
    
    %find region perimeters
    perimlist=bwboundaries(L_final,4,'noholes');
    
    if length(perimlist)~=length(Ctemp)
        warning('Object number mismatch detected during segmentation.')
    end
    
    %prepare output region properties
    area_vec=cat(1,Ctemp.Area);
    theta_vec=cat(1,Ctemp.Orientation)/180*pi;
    C2=cat(1,Ctemp.Centroid);
    Xcent_vec=C2(:,1);
    Ycent_vec=C2(:,2);
    
    %save region properties to data structures, and make sure regions are not perimeter only
    segdata.num_objs=size(C2,1);
    for i=1:segdata.num_objs
            segdata.object(i).Xcent=Xcent_vec(i);
            segdata.object(i).Ycent=Ycent_vec(i);
            segdata.object(i).area=area_vec(i);
            segdata.object(i).theta=theta_vec(i);
            segdata.object(i).Xperim=perimlist{i}(:,2);
            segdata.object(i).Yperim=perimlist{i}(:,1);
        %{  
        %old version 9/22/2016
        if area_vec(i)~=size(perimlist{i},1)
            segdata.object(i).Xcent=Xcent_vec(i);
            segdata.object(i).Ycent=Ycent_vec(i);
            segdata.object(i).area=area_vec(i);
            segdata.object(i).theta=theta_vec(i);
            segdata.object(i).Xperim=perimlist{i}(:,2);
            segdata.object(i).Yperim=perimlist{i}(:,1);
        else
            segdata.object(i).area=0;
        end
    %}
    end
else
    segdata=[];
end

%{
bck1=mat2gray(~(L_final>0).*~G_fill2.*I1);
Ic(:,:,1)=bck1+0.8*abs(G_fill2-(L_final>0));
Ic(:,:,2)=bck1+0.8*(L_final>0);
Ic(:,:,3)=bck1;
imagesc(Ic)
hold on
for j=1:max(L_final(:))
    text(X_cent_vec(j),Y_cent_vec(j),num2str(j),'color',[0 0 1])
end
axis equal tight
box on
title([num2str(N) ' objects'])
%}

%**************************************************************************
% Included functions
%**************************************************************************
%Tristan Ursell
%Relative Noise Transform
%(c) November 2012
%
%Iout=relnoise(Iin,sz,sigma);
%Iout=relnoise(Iin,sz,sigma,'field');
%[Iout,Ivar]=relnoise(Iin,sz,sigma,...);
%[Iout,Ivar,Imean]=relnoise(Iin,sz,sigma,...);
%
%Iin = the input image, of any numerical class.
%
%sz = (3 < sz < min(size(Iin))) is the size of the filter block used to
%calculate means and variances.  This value must be odd.
%
%sigma (sigma > 0) is the weighting parameter that defines the standard
%deviation relative to the filter block's standard deviation around which
%the center pixel will be Gaussian weighted. Setting sigma = 1 weights the
%current pixel using the STD of the current filter block. Lower values
%bring the current pixel closer to the mean, while high values are more
%tolerant of variations.  As sigma -> Inf, Iout = Iin.
%
%The field 'plot' will create an output plot comparing this transform to
%the original image, a Gaussian blur with STD = sz/2, and median filter
%with block size equal to sz.  At first glance, this filter appears similar
%to a median transform, but it does a better job of preserving local
%intensity extrema.  Comparison with the median filter requires the Image
%Processing Toolbox, but the rest of the script does not.
%
%The field 'disk' or 'square' will choose between using a disk or square
%filter block shape, where sz is the disk diameter or square side length.
%The default is square.
%
%The field 'custom' may be followed by a user-defined logical matrix or
%strel, e.g. relnoise(Iin,sz,sigma,'custom',strel('line',30,45)).  In this
%case 'sz' will be unused.
%
%Iout is the transformed output image.
%
%Ivar is the variance of the pixel intensities in the filter block at every
%point in the image -- essentially the spatially varying variance of the
%image.
%
%Imean is the mean smoothed image using the filter block, equivalent to a
%convolution averaging filter with the specified neighborhood.
%
%see also: wiener2  filter2
%
%Example:
%
%Iin=imread('spot_test.tif');
%
%Iout=relnoise(Iin,3,0.5,'square','plot');
% %OR
%Iout=relnoise(Iin,3,0.5,'custom',strel('line',30,0),'plot');
%
%figure;
%subplot(1,2,1)
%imagesc(Iout-double(Iin))
%title('What was removed from the original image.')
%axis equal tight
%box on
%
%subplot(1,2,2)
%imagesc(abs(fftshift(fft(Iout-double(Iin)))))
%title('FFT of difference between original and filtered images.')
%axis equal tight
%box on
%

function [varargout]=relnoise(Iin,sz,sigma,varargin)

%Nans
Inan=isnan(Iin);
if sum(Inan(:))>0
    error('Input matrix contains NaNs.')
end

%convert type
Iin=double(Iin);

%check filter size
if or(sz<1,sz>min(size(Iin)))
    error('The filter size is out of bounds.')
end

if mod(sz,2)~=1
    sz=sz+1;
end

%parse field input
f1=find(strcmp('plot',varargin),1);
f2=find(strcmp('disk',varargin),1);
f3=find(strcmp('custom',varargin),1);

%choose plot option
if ~isempty(f1)
    plotq=1;
else
    plotq=0;
end

%choose filter type
if ~isempty(f3)
    hood_temp=varargin{f3+1};
    if strcmp(class(hood_temp),'strel')
        hood=hood_temp.getnhood;
    else
        hood=hood_temp;
    end
    
    sz=round(mean(size(hood)));
elseif ~isempty(f2)
    %disk filter block
    Xdisk = ones(sz,1)*(-(sz-1)/2:(sz-1)/2);
    Ydisk = (-(sz-1)/2:(sz-1)/2)'*ones(1,sz);
    Zdisk = sqrt(Xdisk.^2 + Ydisk.^2);
    
    hood=zeros(sz,sz);
    hood(Zdisk<=(sz-1)/2)=1;
else
    %square filter block
    hood=ones(sz,sz);
end

%convert hood class
hood=single(hood);

%calcualte means and variances
hood_sz=sum(hood(:));

%perform convolultion and normalization
Imean0=conv2(Iin,hood,'same')/hood_sz;
Inorm=conv2(ones(size(Iin)),hood,'same')/hood_sz;
Imean=Imean0./Inorm;
Ivar=conv2(Iin.^2,hood,'same')./Inorm*1/hood_sz-Imean.^2;
Ivar(Ivar<0)=0;

%compute weight matrix
W=exp(-(Iin-Imean).^2./(2*sigma^2*Ivar));

%correct for zero variance pixels
W(Ivar==0)=0;

%compute output image
if sigma==0
    Iout=Imean;
else
    Iout=Iin.*W+(1-W).*Imean;
end

%handle outputs
if nargout==1
    varargout{1}=Iout;
elseif nargout==2
    varargout{1}=Iout;
    varargout{2}=Ivar;
elseif nargout==3
    varargout{1}=Iout;
    varargout{2}=Ivar;
    varargout{3}=Imean;
elseif nargout==0
else
    error('Incorrect number of output arguments.')
end

%plot comparisons
if plotq==1
    
    figure;
    subplot(2,2,1)
    imagesc(Iin)
    xlabel('X')
    ylabel('Y')
    box on
    axis equal tight
    title('Original Image')
    
    subplot(2,2,2)
    imagesc(Iout)
    xlabel('X')
    ylabel('Y')
    box on
    axis equal tight
    title(['Relative Noise Reduction (this filter), size = ' num2str(sz) ', sigma = ' num2str(sigma)])
    
    subplot(2,2,3)
    %look for image processing toolbox
    boxes=ver;
    gotit=strfind([boxes.Name],'Image Processing Toolbox');
    if ~isempty(gotit)
        Iout2=medfilt2(Iin,[sz,sz],'symmetric');
    else
        disp('Sorry, you do not have the Image Processing Toolbox.')
        disp('The medfilt2 comparison image cannot be generated.')
        disp('Disabling `plot` will stop this message.')
        Iout2=Imean;
    end
    imagesc(Iout2)
    xlabel('X')
    ylabel('Y')
    box on
    axis equal tight
    
    if ~isempty(gotit)
        title(['Median Filter of size ' num2str(sz)])
    else
        title(['Mean Filter of size ' num2str(sz)])
    end
    
    %construct Gaussian filter (without Image Processing Toolbox)
    sz2=round(3/2*sz);
    Xgauss = ones(sz2,1)*(-(sz2-1)/2:(sz2-1)/2);
    Ygauss = (-(sz2-1)/2:(sz2-1)/2)'*ones(1,sz2);
    Zgauss = exp(-(Xgauss.^2+Ygauss.^2)/(2*(sz/2)^2));
    Zgauss = Zgauss/sum(Zgauss(:));
    Iout3=conv2(Iin,Zgauss,'same');
    
    subplot(2,2,4)
    imagesc(Iout3)
    xlabel('X')
    ylabel('Y')
    box on
    axis equal tight
    title(['Gaussian Blur with STD = ' num2str(sz/2)])
    pause(0.1)
end



%Tristan Ursell
%April 2012
%Function to remove false positive regions from segmented bacterial images
%
%Lin image must be bwlabel or segmented BW image.
%Iin image is the grayscale from which Lin was generated.
%hmin sets the ratio of the STD / mean pixel intensities, above which the
%region is a false positive.
%int_rej, together with hmin, sets the ratio of border pixel to interior
%pixel intensities that are the second condition for false positive
%rejection.

function Lout=falsepos(Lin,Iin,int_rej,hmin)

%check intensities of regions between border and interior
L1=bwlabel(Lin);
L1_perim=bwperim(L1);

%find the list of regions that appear in both interior and perimeter forms
%this is to handle the very rare case that a region only has perimeter
%pixels
L1_list=L1-L1_perim.*L1;
v1=unique(L1(L1>0));
v2=unique(L1_list(L1_list>0));
setdiff1=setdiff(v1,v2);

for i=1:length(setdiff1)
    L1(L1==setdiff1(i))=0;
end

%reform inputs
L1=bwlabel(L1);
L1_perim=bwperim(L1);
Lout=L1;

if max(L1(:))==0
    return
end

%find statistics for segmented regions
temp1=regionprops(~L1_perim.*L1,Iin,'MeanIntensity','PixelValues');
temp2=regionprops(L1_perim.*L1,Iin,'MeanIntensity');

interiormean=[temp1.MeanIntensity]; %segmented interior pixel means
perimmean=[temp2.MeanIntensity]; %segemented border pixel means

interiorstd=zeros(1,length(temp1)); %segmented interior pixel STDs
for j=1:max(L1(:))
    interiorstd(j)=std(temp1(j).PixelValues);
end

%calculate border / intertior ratios (should be greater than 1, unless peripheral)
brd_ratio=perimmean./interiormean;

%nois / mean ratio
noise_ratio=interiorstd./interiormean;

%remove false positive regions based on border / interior ratios
cond1=brd_ratio<=(1+int_rej*hmin); %if the mean of the region perimeter pixels is inappropriate compared to the mean of the interior pixels
cond2=noise_ratio<hmin/2; %if the STD of the interior pixel regions is too small
false_pos=find(or(cond1,cond2));
for j=1:length(false_pos)
    Lout(L1==false_pos(j))=0;
end

%turn to logical, and remove spurs
Lout=bwmorph(Lout>0,'spur',Inf);










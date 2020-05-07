function varargout = view_contours(varargin)
% VIEW_CONTOURS MATLAB code for view_contours.fig
%      VIEW_CONTOURS, by itself, creates a new VIEW_CONTOURS or raises the existing
%      singleton*.
%
%      H = VIEW_CONTOURS returns the handle to a new VIEW_CONTOURS or the handle to
%      the existing singleton*.
%
%      VIEW_CONTOURS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_CONTOURS.M with the given input arguments.
%
%      VIEW_CONTOURS('Property','Value',...) creates a new VIEW_CONTOURS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_contours_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_contours_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_contours

% Last Modified by GUIDE v2.5 22-Jun-2015 23:31:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @view_contours_OpeningFcn, ...
    'gui_OutputFcn',  @view_contours_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before view_contours is made visible.
function view_contours_OpeningFcn(hObject, eventdata, handles, varargin)
global Nim
global im_name
global cnt_name
global frame
global Ncell
global c_mat
global make_eps

Nim=1;  %number of images in the stack
im_name=[]; %full path and file name of image stack
cnt_name=[]; %full path and file name of contour *.mat file
frame=[]; %'frame' data structure loaded from mat file
Ncell=0; %number of labeled cells in mat file
c_mat=[]; %distinct colors for contours
make_eps=0;

set(handles.checkbox_acrosstime,'Visible','off')

handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = view_contours_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_viewcontours.
function pushbutton_viewcontours_Callback(hObject, eventdata, handles)
global Nim
global im_name
global cnt_name
global frame
global Ncell
global c_mat

%get flr stack name
[file1,path1]=uigetfile({'*.tif;*.tiff','TIF Image Files'},'Select Signal Stack');
if file1==0
    return
end
cd(path1)

[file2,path2]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
if file2==0
    return
end
cd(path2)

%construct names
im_name=fullfile(path1,file1);
cnt_name=fullfile(path2,file2);

set(handles.text_contourfile,'String',' ')
set(handles.text_imagefile,'String','Loading mat-file ...')
drawnow 

load(fullfile(path2,file2),'frame','Ncell')

if ~isfield(frame(1).object,'Xcont')
    set(handles.text_imagefile,'String',['No contour data found.'])
    pause(2)
end

set(handles.text_imagefile,'String',[file1])
set(handles.text_contourfile,'String',[file2])

%get number of images
Nim=length(frame);
set(handles.text_frameend,'String',[' of ' num2str(Nim)])

%choose contour colors (if requested)
c_temp=jet(Ncell);
[~,rand_vec]=sort(rand(Ncell,1));
c_mat=c_temp(rand_vec,:);
c_mat(c_mat<=0)=0;
c_mat(c_mat>=1)=1;


% --- Executes on button press in pushbutton_process.
function pushbutton_process_Callback(hObject, eventdata, handles)
global Nim
global im_name
global cnt_name
global frame
global c_mat

if isempty(im_name)
    warndlg('First load image and contour data.','Load data first');
    return
end

%pole markersize on plot
mk_sz=15;

%set text number color
num_color=[0,1,1];

%load original image
I0=zeros(size(imread(im_name,1)));

%get contour shift values
dX=str2num(get(handles.edit_Ch1_X,'String'));
dY=str2num(get(handles.edit_Ch1_Y,'String'));

%remove previous  file if saving TIF and not using arrow keys
[path1,fname1,~]=fileparts(cnt_name);     
%tif_name=fullfile(path1,[date '_' fname1 '.tif']);
tif_name=fullfile(path1,['Tracking_' fname1 '.tif']);

if and(~get(handles.checkbox_keypress,'Value'),get(handles.checkbox_savetiff,'Value'))
    delete(tif_name)
end

%frame number
f1=figure('units','normalized','outerposition',[0 0 1 1]);
a1 = axes('Parent',f1);

i=0;
while i<Nim
    i=i+1;
    
    %load original image
    Im0=mat2gray(imread(im_name,i));
    
    %background removal
    if exist('back_mean','var')
        Im0=Im0-back_mean(i);
    end
    
    %plotting
    if get(handles.checkbox_keypress,'Value')
        pause(str2num(get(handles.edit_framedelay,'String')))
    end
    
    figure(f1);
    hold off
    imagesc(Im0)
    colormap(gray)
    hold on
    xlabel('X')
    ylabel('Y')
    axis equal tight
    title(['Frame: ' num2str(i) ', Cells: ' num2str(frame(i).num_objs)])
    
    
    tic
    nocont=0;
    for j=1:frame(i).num_objs
        if and(isfield(frame(i).object(j),'cellID'),isfield(frame(i).object(j),'Xcont'))
            if and(~isempty(frame(i).object(j).Xcont),~isempty(frame(i).object(j).cellID))
                %get current contour points
                X=frame(i).object(j).Xcont;
                Y=frame(i).object(j).Ycont;
                
                %perform contour shift
                X=X+dX;
                Y=Y+dY;
                
                %choose contour color
                if get(handles.checkbox_unqcolor,'Value')
                    c_vec=c_mat(frame(i).object(j).cellID,:);
                else %non-unique color
                    if get(handles.checkbox_morphoclass,'Value')
                        c_vec=[0.1,0.1,0.1];
                    else
                        c_vec=[0,0.9,0];
                    end
                end
                
                %get mesh if available
                if and(get(handles.checkbox_showmesh,'Value'),...
                        and(~isfield(frame(i).object(j),'pill_mesh'),...
                        ~isfield(frame(i).object(j),'mesh')))
                    disp(['Mesh plotting requested, but no mesh data found.'])
                end
                
                if get(handles.checkbox_showmesh,'Value')
                    if isfield(frame(i).object(j),'pill_mesh')
                        mesh_mat=frame(i).object(j).pill_mesh;
                        
                        for k=1:size(mesh_mat,1)
                            line([mesh_mat(k,1),mesh_mat(k,3)],[mesh_mat(k,2),mesh_mat(k,4)],'color',c_vec,'Parent',a1)
                        end
                        
                    elseif isfield(frame(i).object(j),'mesh')
                        for k=1:length(frame(i).object(j).mesh)
                            line(frame(i).object(j).mesh(k).X_pts,frame(i).object(j).mesh(k).Y_pts,'color',c_vec,'Parent',a1)
                        end
                    end
                end
                
                %plot the contour
                line(X,Y,'LineWidth',1,'Color',c_vec,'Parent',a1);
                
                %plot bulges, constrictions and bends
                if and(get(handles.checkbox_morphoclass,'Value'),isfield(frame(i).object(j),'pill_mesh'))
                    %get the match contour curvatures
                    kappa1=frame(i).object(j).side1_kappa;
                    kappa2=frame(i).object(j).side2_kappa;
                    
                    mesh_mat=frame(i).object(j).pill_mesh;
                    
                    for k=1:length(kappa1)
                        xtemp=[mesh_mat(k,1),mesh_mat(k,3)];
                        ytemp=[mesh_mat(k,2),mesh_mat(k,4)];
                        
                        if and(kappa1(k)>=0,kappa2(k)>=0) %bulge
                            line(xtemp,ytemp,'color',[1 0.2 0],'Parent',a1)
                        end
                        if and(kappa1(k)<=0,kappa2(k)<=0) %constriction
                            line(xtemp,ytemp,'color',[0 0.3 1],'Parent',a1)
                        end
                        if or(and(kappa1(k)>0,kappa2(k)<0),and(kappa1(k)<0,kappa2(k)>0)) %bend
                            line(xtemp,ytemp,'color',[0 0.9 0],'Parent',a1)
                        end
                        if and(kappa1(k)==0,kappa2(k)==0) %straight
                            line(xtemp,ytemp,'color',[0 0 0],'Parent',a1)
                        end
                    end
                end
                
                if get(handles.checkbox_polemark,'Value')
                    %get poles
                    p1=frame(i).object(j).pole1;
                    p2=frame(i).object(j).pole2;
                    
                    plot(X(1),Y(1),'.','color',[0 1 0],'markersize',mk_sz,'Parent',a1)
                    plot(X(end),Y(end),'.','color',[1 0 0],'markersize',mk_sz,'Parent',a1)
                    
                    if p1==1
                        plot(X(p2),Y(p2),'rx','markersize',7,'Parent',a1)
                    else
                        plot(X(p1),Y(p1),'rx','markersize',7,'Parent',a1)
                    end
                end
                
                if get(handles.checkbox_numbers,'Value')
                    if and(get(handles.checkbox_unqcolor,'Value'),get(handles.checkbox_showmesh,'Value'))
                        text(mean(X),mean(Y),num2str(frame(i).object(j).cellID),'color',([1,1,1]-c_vec).^2,'Parent',a1)
                    else
                        text(mean(X),mean(Y),num2str(frame(i).object(j).cellID),'color',num_color,'Parent',a1)
                    end
                end
            end
        end
        
        if ~isfield(frame(i).object(j),'Xcont')
            if ~nocont
                nocont=1;
                mask_seg=zeros(size(Im0));
                im_sz=size(Im0);
            end
            
            %convert perimeter pixels into 1D array and turn them white 
            perim_inds=sub2ind(im_sz,frame(i).object(j).Yperim,frame(i).object(j).Xperim);
            mask_seg(perim_inds)=1;
        end
    end
    
    if exist('mask_seg','var')
        mask_seg=imfill(mask_seg,'holes');
        
        T(:,:,1)=Im0.*~mask_seg;
        T(:,:,2)=Im0;
        T(:,:,3)=Im0.*~mask_seg;
        
        hold off
        imagesc(T)
        hold on
        xlabel('X')
        ylabel('Y')
        axis equal tight
        title(['Frame: ' num2str(i) ', Cells: ' num2str(frame(i).num_objs)])
        
        for j=1:frame(i).num_objs
            if get(handles.checkbox_numbers,'Value')
                text(mean(frame(i).object(j).Xperim),mean(frame(i).object(j).Yperim),num2str(frame(i).object(j).cellID),'color',[1 1 1]-num_color,'Parent',a1)
            end
        end
        clearvars mask_seg
    end
    drawnow
    dt1=toc;
    
    %handle arrow key presses
    if get(handles.checkbox_keypress,'Value')
        title(['Frame: ' num2str(i) ', Cells: ' num2str(frame(i).num_objs) ', use left and right arrow keys to progress, esc to exit'])
        %set(f1,'KeyPressFcn','key_pressed')            
        set(f1,'KeyPressFcn') %Fixes the problem - added by Andres Florez on 04/16/2020 
        try
            waitforbuttonpress
        catch
            return
        end
        key=get(f1,'CurrentKey');
        
        %[~,~,key]=ginput(1);
        
        while all([~strcmp(key,'escape'),~strcmp(key,'leftarrow'),~strcmp(key,'rightarrow')])
            %[~,~,key]=ginput(1);
            %set(f1,'KeyPressFcn','key_pressed')
            set(f1,'KeyPressFcn') %Fixes the problem - added by Andres Florez on 04/16/2020 
            waitforbuttonpress
            key=get(f1,'CurrentKey');
        end
        
        %forward arrow
        if strcmp(key,'rightarrow')
        end
        
        %back arrow
        if and(i>1,strcmp(key,'leftarrow'))
            i=i-2;
        elseif and(i==1,strcmp(key,'leftarrow'))
            i=i-1;
        end
        
        %esc key, stop everything
        if strcmp(key,'escape')
            close(f1)
            return
        end
    else
        pause((str2num(get(handles.edit_framedelay,'String'))-dt1).*((num2str(get(handles.edit_framedelay,'String'))-dt1)>0))
    end

%end

%Save tif stack function (below) has to be included in the while loop, otherwise it doesn't
%save the whole tif stack but only the last frame. -Andres

if and(get(handles.checkbox_savetiff,'Value'),~get(handles.checkbox_keypress,'Value'))

       % added by Andres Florez on 04/16/2020
       print(f1,'-dtiff','-r60','contour_temp.tif');%if you change resolution, change cropping box
       %print(f1,'-dtiff','-r1400','contour_temp.tif')
       Itemp=imread('contour_temp.tif');
       %to make the crop box, run [J, rect] = imcrop(Itemp) select a
       %square, hit enter and take the rect vector and paste it here
       %Itempsmall=imcrop(Itemp,[745.51 113.51 2114.98 1585.98]);
       Itempsmall=imcrop(Itemp,[287.51 44.51 938.98 708.98]);
       image_save(Itempsmall,tif_name)
       delete contour_temp.tif
       
       %Old code:
       %export_fig(f1,tif_name,'-tiff' ,'-r100','-append') % this solution is faster - Andres
       %export_fig(f1,'contour_temp.tif','-tiff' ,'-r100')
       %Itempmerge= cat(3, Itemp(:,:,1),Itemp(:,:,2),Itemp(:,:,3));
       %image_save(Itemp,tif_name)
       %export_fig(Itempmerge,tif_name,'-tiff' ,'-r100')
end

end %Add the end here to fix the bug - added by Andres Florez on 04/16/2020

% --- Executes on button press in pushbutton_viewframe.
function pushbutton_viewframe_Callback(hObject, eventdata, handles)
global im_name
global frame
global c_mat
global make_eps

if isempty(im_name)
    warndlg('First load image and contour data.','Load data first');
    return
end

%pole markersize on plot
mk_sz=15;

%set text number color
num_color=[0,1,1];

%get contour shift values
dX=str2num(get(handles.edit_Ch1_X,'String'));
dY=str2num(get(handles.edit_Ch1_Y,'String'));

%frame number
f1=figure('units','normalized','outerposition',[0 0 1 1]);
i=str2num(get(handles.edit_viewframe,'String'));

%load original image
Im0=mat2gray(imread(im_name,i));

%plotting
figure(f1);
a1 = axes('Parent',f1);
hold off
imagesc(Im0)
hold on
colormap(gray)
title(['Frame: ' num2str(i) ', Cells: ' num2str(frame(i).num_objs)])
xlabel('X')
ylabel('Y')
axis equal tight

nocont=0;
for j=1:frame(i).num_objs
    if and(isfield(frame(i).object(j),'cellID'),isfield(frame(i).object(j),'Xcont'))
        if and(~isempty(frame(i).object(j).Xcont),~isempty(frame(i).object(j).cellID))
            %get current contour points
            X=frame(i).object(j).Xcont;
            Y=frame(i).object(j).Ycont;
            
            %perform contour shift
            X=X+dX;
            Y=Y+dY;
            
            %choose contour color
            if get(handles.checkbox_unqcolor,'Value')
                c_vec=c_mat(frame(i).object(j).cellID,:);
            else %non-unique color
                if get(handles.checkbox_morphoclass,'Value')
                    c_vec=[0.1,0.1,0.1];
                else
                    c_vec=[0,0.9,0];
                end
            end
            
            %get mesh if available
            if and(get(handles.checkbox_showmesh,'Value'),...
                    and(~isfield(frame(i).object(j),'pill_mesh'),...
                    ~isfield(frame(i).object(j),'mesh')))
                display(['Mesh plotting requested, but no mesh data found.'])
            end
            
            if get(handles.checkbox_showmesh,'Value')
                if and(isfield(frame(i).object(j),'pill_mesh'),~get(handles.checkbox_morphoclass,'Value'))
                    mesh_mat=frame(i).object(j).pill_mesh;
                    
                    for k=1:size(mesh_mat,1)
                        line([mesh_mat(k,1),mesh_mat(k,3)],[mesh_mat(k,2),mesh_mat(k,4)],'color',c_vec,'Parent',a1)
                    end
                    
                    line([frame(i).object(j).centerline(:,1)],[frame(i).object(j).centerline(:,2)],'color',c_vec,'Parent',a1)
                    
                elseif isfield(frame(i).object(j),'mesh')
                    for k=1:length(frame(i).object(j).mesh)
                        line(frame(i).object(j).mesh(k).X_pts,frame(i).object(j).mesh(k).Y_pts,'color',c_vec,'Parent',a1)
                    end
                end
            end
            
            %plot the contour
            line(X,Y,'LineWidth',1,'Color',c_vec,'Parent',a1);
            
            %plot bulges, constrictions and bends
            if and(get(handles.checkbox_morphoclass,'Value'),isfield(frame(i).object(j),'pill_mesh'))
                %get the match contour curvatures
                kappa1=frame(i).object(j).side1_kappa;
                kappa2=frame(i).object(j).side2_kappa;
                
                mesh_mat=frame(i).object(j).pill_mesh;
                
                for k=1:length(kappa1)
                    xtemp=[mesh_mat(k,1),mesh_mat(k,3)];
                    ytemp=[mesh_mat(k,2),mesh_mat(k,4)];
                    
                    if and(kappa1(k)>=0,kappa2(k)>=0) %bulge
                        line(xtemp,ytemp,'color',[1 0.2 0],'Parent',a1)
                    end
                    if and(kappa1(k)<=0,kappa2(k)<=0) %constriction
                        line(xtemp,ytemp,'color',[0 0.3 1],'Parent',a1)
                    end
                    if or(and(kappa1(k)>0,kappa2(k)<0),and(kappa1(k)<0,kappa2(k)>0)) %bend
                        line(xtemp,ytemp,'color',[0 0.9 0],'Parent',a1)
                    end
                    if and(kappa1(k)==0,kappa2(k)==0) %straight
                        line(xtemp,ytemp,'color',[0 0 0],'Parent',a1)
                    end
                end
            end
            
            if get(handles.checkbox_polemark,'Value')
                %get poles
                p1=frame(i).object(j).pole1;
                p2=frame(i).object(j).pole2;
                
                plot(X(1),Y(1),'.','color',[0 1 0],'markersize',mk_sz,'Parent',a1)
                plot(X(end),Y(end),'.','color',[1 0 0],'markersize',mk_sz,'Parent',a1)
                
                if p1==1
                    plot(X(p2),Y(p2),'rx','markersize',7,'Parent',a1)
                else
                    plot(X(p1),Y(p1),'rx','markersize',7,'Parent',a1)
                end
            end
            
            if get(handles.checkbox_numbers,'Value')
                if and(get(handles.checkbox_unqcolor,'Value'),get(handles.checkbox_showmesh,'Value'))
                    text(mean(X),mean(Y),num2str(frame(i).object(j).cellID),'color',([1,1,1]-c_vec).^2,'Parent',a1)
                else
                    text(mean(X),mean(Y),num2str(frame(i).object(j).cellID),'color',num_color,'Parent',a1)
                end
            end
        end
    end
    if ~isfield(frame(i).object(j),'Xcont')
        if ~nocont
            nocont=1;
            mask_seg=zeros(size(Im0));
            im_sz=size(Im0);
        end
        
        %convert perimeter pixels into 1D array and turn them white
        perim_inds=sub2ind(im_sz,frame(i).object(j).Yperim,frame(i).object(j).Xperim);
        mask_seg(perim_inds)=1;
    end
end

if exist('mask_seg','var')
    mask_seg=imfill(mask_seg,'holes');
    
    T(:,:,1)=Im0.*~mask_seg;
    T(:,:,2)=Im0;
    T(:,:,3)=Im0.*~mask_seg;
    
    hold off
    imagesc(T)
    hold on
    xlabel('X')
    ylabel('Y')
    axis equal tight
    title(['Frame: ' num2str(i) ', Cells: ' num2str(frame(i).num_objs)])
    
    for j=1:frame(i).num_objs
        if get(handles.checkbox_numbers,'Value')
            text(mean(frame(i).object(j).Xperim),mean(frame(i).object(j).Yperim),num2str(frame(i).object(j).cellID),'color',[1 1 1]-num_color,'Parent',a1)
        end
    end
    clearvars mask_seg
end
drawnow

if make_eps
    make_eps=0;
    [path1,file1,~]=fileparts(im_name);
    print(f1,'-depsc2',[fullfile(path1,file1) '_' date '_frame_' num2str(i) '.eps'])
end

% --- Executes on button press in checkbox_filter.
function checkbox_polemark_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_savetiff.
function checkbox_savetiff_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_numbers.
function checkbox_numbers_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_showcolorcurve.
function checkbox_showcolorcurve_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_unqcolor.
function checkbox_unqcolor_Callback(hObject, eventdata, handles)
set(handles.checkbox_morphoclass,'Value',0)

% --- Executes on button press in checkbox_morphoclass.
function checkbox_morphoclass_Callback(hObject, eventdata, handles)
set(handles.checkbox_unqcolor,'Value',0)
set(handles.checkbox_showmesh,'Value',0)

% --- Executes on button press in checkbox_showmesh.
function checkbox_showmesh_Callback(hObject, eventdata, handles)
set(handles.checkbox_morphoclass,'Value',0)

function edit_Ch1_X_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_Ch1_X_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Ch1_Y_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_Ch1_Y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_viewframe_Callback(hObject, eventdata, handles)
global Nim
temp1=str2num(get(handles.edit_viewframe,'String'));
if temp1<1
    temp1=1;
elseif temp1>Nim
    temp1=Nim;
end
set(handles.edit_viewframe,'String',num2str(temp1))
% --- Executes during object creation, after setting all properties.
function edit_viewframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_contourplots.
function pushbutton_contourplots_Callback(hObject, eventdata, handles)
global frame
global Nim

%check if data has been loaded
if isempty(frame)
    [file2,path2]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
    if file2==0
        return
    end
    cd(path2)
    
    %construct names
    set(handles.text_contourfile,'String',' ')
    set(handles.text_imagefile,'String','Loading mat-file ...')
    drawnow
    
    %load data
    load(fullfile(path2,file2),'frame','Ncell');
    set(handles.text_imagefile,'String',' ')
    
    set(handles.text_contourfile,'String',[file2])
    
    %get number of images
    Nim=length(frame);
    set(handles.text_frameend,'String',[' of ' num2str(Nim)])
elseif ~isempty(frame)
    if ~isfield(frame(1).object(1),'pill_mesh')
        [file2,path2]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
        if file2==0
            return
        end
        cd(path2)
        
        %construct names
        set(handles.text_contourfile,'String',' ')
        set(handles.text_imagefile,'String','Loading mat-file ...')
        drawnow
        
        %load data
        load(fullfile(path2,file2),'frame','Ncell');
        set(handles.text_imagefile,'String',' ')
        
        set(handles.text_contourfile,'String',[file2])
        
        %get number of images
        Nim=length(frame);
        set(handles.text_frameend,'String',[' of ' num2str(Nim)])
    end
end

%check for correct file type
if ~isfield(frame(1).object(1),'pill_mesh')
    warndlg('The file you selected does not contain usable data.')
end

%get all input parameters
dt=str2num(get(handles.edit_frametime,'String'));
calb=str2num(get(handles.edit_plotcalb,'String'));
nbin=str2num(get(handles.edit_binnum,'String'));
curv_cut=str2num(get(handles.edit_curvcut,'String'));

dt_vec=linspace(0,(Nim-1)*dt,Nim);
    
%get checkbox values
chbx1=get(handles.checkbox_widthhist,'Value');
chbx2=get(handles.checkbox_lengthhist,'Value');
chbx3=get(handles.checkbox_meanwidth,'Value');
chbx4=get(handles.checkbox_meanlength,'Value');
chbx_std=get(handles.checkbox_wstd,'Value');
chbx_prox=get(handles.checkbox_exproximal,'Value');

%width and length distributions
if or(chbx1,chbx2)
    %without time
    nums=[];
    
    %handle length
    if chbx2
        for i=1:length(frame)
            %raw data
            data1=[frame(i).object.length]*calb;
            
            %test if cells are proximal
            for j=1:frame(i).num_objs
                if chbx_prox
                    condprox(j)=isempty(frame(i).object(j).proxID);
                else
                    condprox(j)=1;
                end
            end
            
            %aggregate data
            nums=[nums,data1(condprox>0)];
        end
        
        if max(nums)>3000
            x_label='Length (um)';
            nums=nums/1000;
        else
            x_label='Length (nm)';
        end
        
        f1=figure;
        hist(nums,nbin)
        xlabel(x_label)
        ylabel('Frequency')
        title(['Total number of objects:  ' num2str(length(nums))])
        box on
    end
    
    %handle width
    if chbx1
        for i=1:length(frame)
            
            %calcuate width
            w1=zeros(1,frame(i).num_objs);
            for j=1:frame(i).num_objs
                if or(~chbx_prox,isempty(frame(i).object(j).proxID))
                    %get base widths and curvatures
                    w0=[frame(i).object(j).width]*calb;
                    kappa1=frame(i).object(j).side1_kappa*1000/calb;
                    kappa2=frame(i).object(j).side2_kappa*1000/calb;
                    
                    %remove endpoints
                    w0=w0(2:end-1);
                    kappa1=kappa1(2:end-1);
                    kappa2=kappa2(2:end-1);
                    
                    %filter based on curvature
                    if ~isempty(curv_cut)
                        cond1=and(kappa1<curv_cut,kappa2<curv_cut);
                    else
                        cond1=logical(ones(size(kappa1)));
                    end
                    
                    %filtered, average width for each object
                    w1(j)=mean(w0(cond1));
                end
            end
            
            %aggregate data
            nums=[nums,w1(w1>0)];
        end
        
        if max(nums)>3000
            x_label='Width (um)';
            nums=nums/1000;
        else
            x_label='Width (nm)';
        end
        
        f1=figure;
        hist(nums,nbin)
        xlabel(x_label)
        ylabel('Frequency')
        title(['Total number of objects:  ' num2str(length(nums))])
        box on
    end
end

%mean width and length plots
if or(chbx3,chbx4)
    x_label='Time (s)';
    
    if max(dt_vec)>300
        dt_vec=dt_vec/60;
        x_label='Time (min)';
    end
    
    for i=1:length(frame)
        %handle length
        if chbx4
            %raw data
            data1=[frame(i).object.length];
            
            %test if cells are proximal
            for j=1:frame(i).num_objs
                if chbx_prox
                    condprox(j)=isempty(frame(i).object(j).proxID);
                else
                    condprox(j)=1;
                end
            end
            
            meanL(i)=mean(data1(condprox>0))*calb;
            stdL(i)=std(data1(condprox>0))*calb;
        end
        
        %handle width
        if chbx3
            w1=zeros(1,frame(i).num_objs);
            for j=1:frame(i).num_objs
                if or(~chbx_prox,isempty(frame(i).object(j).proxID))
                    %get base widths and curvatures
                    w0=[frame(i).object(j).width]*calb;
                    kappa1=frame(i).object(j).side1_kappa*1000/calb;
                    kappa2=frame(i).object(j).side2_kappa*1000/calb;
                    
                    %remove endpoints
                    w0=w0(2:end-1);
                    kappa1=kappa1(2:end-1);
                    kappa2=kappa2(2:end-1);
                    
                    %filter based on curvature
                    if ~isempty(curv_cut)
                        cond1=and(kappa1<curv_cut,kappa2<curv_cut);
                    else
                        cond1=logical(ones(size(kappa1)));
                    end
                    
                    %filtered, average width for each object
                    w1(j)=mean(w0(cond1));
                end
            end
            
            %mean of whole cell widths
            meanW(i)=mean(w1(w1>0));
            stdW(i)=std(w1(w1>0));
        end
    end
    
    %plot mean length
    if chbx4
        if max(meanL)>3000
            y_label='Length (um)';
            meanL=meanL/1000;
            stdL=stdL/1000;
        else
            y_label='Length (nm)';
        end
        
        f1=figure;
        hold on
        if chbx_std
            %patch([dt_vec,flip(dt_vec)],[meanL+stdL,flip(meanL-stdL)],[0.8 0.8 0.8],'EdgeColor','none')
            fill([dt_vec,flip(dt_vec)],[meanL+stdL,flip(meanL-stdL)],[0.8 0.8 0.8],'EdgeColor','none')
        end
        plot(dt_vec,meanL,'b-','linewidth',1.5)
        xlabel(x_label)
        ylabel(y_label)
        dtt=(max(dt_vec)-min(dt_vec))/100;
        xlim([min(dt_vec)-dtt,max(dt_vec)+dtt])
        box on
    end
    
    %plot mean width
    if chbx3
        if max(meanW)>3000
            y_label='Width (um)';
            meanW=meanW/1000;
            stdW=stdW/1000;
        else
            y_label='Width (nm)';
        end
        
        f1=figure;
        hold on
        if chbx_std
            %patch([dt_vec,flip(dt_vec)],[meanL+stdL,flip(meanL-stdL)],[0.8 0.8 0.8],'EdgeColor','none')
            fill([dt_vec,flip(dt_vec)],[meanW+stdW,flip(meanW-stdW)],[0.8 0.8 0.8],'EdgeColor','none')
        end
        plot(dt_vec,meanW,'b-','linewidth',1.5)
        xlabel(x_label)
        ylabel(y_label)
        dtt=(max(dt_vec)-min(dt_vec))/100;
        xlim([min(dt_vec)-dtt,max(dt_vec)+dtt])
        box on
    end
end


%*************************************************************
%*************************************************************
%*************************************************************


% --- Executes on button press in pushbutton_frame_eps.
function pushbutton_frame_eps_Callback(hObject, eventdata, handles)
global make_eps
make_eps=1;
pushbutton_viewframe_Callback(hObject, eventdata, handles)

function edit_framedelay_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_framedelay,'String'))<0
    set(handles.edit_framedelay,'String','0')
end
% --- Executes during object creation, after setting all properties.
function edit_framedelay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_keypress.
function checkbox_keypress_Callback(hObject, eventdata, handles)
if get(handles.checkbox_keypress,'Value')
    set(handles.text11,'enable','off')
    set(handles.edit_framedelay,'enable','off')
    set(handles.checkbox_savetiff,'enable','off')
else
    set(handles.text11,'enable','on')
    set(handles.edit_framedelay,'enable','on')
    set(handles.checkbox_savetiff,'enable','on')
end


% --- Executes on button press in checkbox_widthhist.
function checkbox_widthhist_Callback(hObject, eventdata, handles)
set(handles.checkbox_widthhist,'Value',1)
set(handles.checkbox_lengthhist,'Value',0)
set(handles.checkbox_meanwidth,'Value',0)
set(handles.checkbox_meanlength,'Value',0)
set(handles.text21,'enable','on')
set(handles.edit_binnum,'enable','on')
set(handles.checkbox_wstd,'enable','off')
set(handles.checkbox_acrosstime,'enable','on')
set(handles.text22,'enable','on')
set(handles.edit_curvcut,'enable','on')

if get(handles.checkbox_acrosstime,'Value')
    set(handles.text20,'enable','on')
    set(handles.edit_frametime,'enable','on')
else
    set(handles.text20,'enable','off')
    set(handles.edit_frametime,'enable','off')
end

% --- Executes on button press in checkbox_lengthhist.
function checkbox_lengthhist_Callback(hObject, eventdata, handles)
set(handles.checkbox_lengthhist,'Value',1)
set(handles.checkbox_widthhist,'Value',0)
set(handles.checkbox_meanwidth,'Value',0)
set(handles.checkbox_meanlength,'Value',0)
set(handles.text21,'enable','on')
set(handles.edit_binnum,'enable','on')
set(handles.checkbox_wstd,'enable','off')
set(handles.checkbox_acrosstime,'enable','on')
set(handles.text22,'enable','off')
set(handles.edit_curvcut,'enable','off')

if get(handles.checkbox_acrosstime,'Value')
    set(handles.text20,'enable','on')
    set(handles.edit_frametime,'enable','on')
else
    set(handles.text20,'enable','off')
    set(handles.edit_frametime,'enable','off')
end

% --- Executes on button press in checkbox_acrosstime.
function checkbox_acrosstime_Callback(hObject, eventdata, handles)
if or(get(handles.checkbox_lengthhist,'Value'),get(handles.checkbox_widthhist,'Value'))
    if get(handles.checkbox_acrosstime,'Value')
        set(handles.text20,'enable','on')
        set(handles.edit_frametime,'enable','on')
    else
        set(handles.text20,'enable','off')
        set(handles.edit_frametime,'enable','off')
    end
end


% --- Executes on button press in checkbox_meanlength.
function checkbox_meanlength_Callback(hObject, eventdata, handles)
set(handles.checkbox_meanlength,'Value',1)
set(handles.checkbox_lengthhist,'Value',0)
set(handles.checkbox_meanwidth,'Value',0)
set(handles.checkbox_widthhist,'Value',0)
set(handles.text20,'enable','on')
set(handles.edit_frametime,'enable','on')
set(handles.text21,'enable','off')
set(handles.edit_binnum,'enable','off')
set(handles.checkbox_wstd,'enable','on')
set(handles.checkbox_acrosstime,'enable','off')
set(handles.text22,'enable','off')
set(handles.edit_curvcut,'enable','off')


% --- Executes on button press in checkbox_meanwidth.
function checkbox_meanwidth_Callback(hObject, eventdata, handles)
set(handles.checkbox_meanwidth,'Value',1)
set(handles.checkbox_lengthhist,'Value',0)
set(handles.checkbox_widthhist,'Value',0)
set(handles.checkbox_meanlength,'Value',0)
set(handles.text20,'enable','on')
set(handles.edit_frametime,'enable','on')
set(handles.text21,'enable','off')
set(handles.edit_binnum,'enable','off')
set(handles.checkbox_wstd,'enable','on')
set(handles.checkbox_acrosstime,'enable','off')
set(handles.text22,'enable','on')
set(handles.edit_curvcut,'enable','on')


% --- Executes on button press in checkbox_wstd.
function checkbox_wstd_Callback(hObject, eventdata, handles)


function edit_frametime_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_frametime,'String')) <0
    set(handles.edit_frametime,'String',num2str(1))
end
% --- Executes during object creation, after setting all properties.
function edit_frametime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_binnum_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_binnum,'String')) <1
    set(handles.edit_binnum,'String',num2str(1))
end
% --- Executes during object creation, after setting all properties.
function edit_binnum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_plotcalb_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_plotcalb,'String')) <0
    set(handles.edit_plotcalb,'String',num2str(1))
end
% --- Executes during object creation, after setting all properties.
function edit_plotcalb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_curvcut_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_curvcut_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_exproximal.
function checkbox_exproximal_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_clearbutton.
function pushbutton_clearbutton_Callback(hObject, eventdata, handles)
global frame
global Nim
clearvars -global frame Nim
set(handles.text_contourfile,'String','Data cleared.')
set(handles.text_imagefile,'String',' ')

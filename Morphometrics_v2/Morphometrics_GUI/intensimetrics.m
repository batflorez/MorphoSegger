function varargout = intensimetrics(varargin)
% INTENSIMETRICS MATLAB code for intensimetrics.fig
%      INTENSIMETRICS, by itself, creates a new INTENSIMETRICS or raises the existing
%      singleton*.
%
%      H = INTENSIMETRICS returns the handle to a new INTENSIMETRICS or the handle to
%      the existing singleton*.
%
%      INTENSIMETRICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTENSIMETRICS.M with the given input arguments.
%
%      INTENSIMETRICS('Property','Value',...) creates a new INTENSIMETRICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before intensimetrics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to intensimetrics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help intensimetrics

% Last Modified by GUIDE v2.5 16-Jan-2016 20:32:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @intensimetrics_OpeningFcn, ...
    'gui_OutputFcn',  @intensimetrics_OutputFcn, ...
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


% --- Executes just before intensimetrics is made visible.
function intensimetrics_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for intensimetrics
handles.output = hObject;

global mat_file
global im_file

%initialize
mat_file=[];
im_file=[];

%default values
set(handles.checkbox_allcells,'Value',1)
set(handles.edit_cellnumber,'Enable','off')
set(handles.checkbox_ignoreproximal,'Value',1)
set(handles.checkbox_background,'Value',1)
set(handles.checkbox_contour,'Value',1)
set(handles.checkbox_align,'Value',1)
set(handles.checkbox_smooth,'Value',1)
set(handles.checkbox_plotq,'Value',1)
set(handles.checkbox_noise,'Value',1)

%defaults for contour intensity profile
global dl
global res
global offset

dl=5;  %kymograph project length (pixels)
res=10;  %profile vector resolution
offset=0.7;  %offset percentage inside

global rel_sz
global rel_sigma

rel_sz=3; %size of the relative noise filter
rel_sigma=1; %sigma parameter for rel-noise filter

set(handles.checkbox_align,'Visible','off') %look in get_kymo_data_from_file_DF.m

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = intensimetrics_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_matfile.
function pushbutton_matfile_Callback(hObject, eventdata, handles)
global mat_file

[file2,path2]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
if file2==0
    return
end
cd(path2)

mat_file=fullfile(path2,file2);

set(handles.text_status,'String',['loaded: ' file2])


% --- Executes on button press in pushbutton_imagestack.
function pushbutton_imagestack_Callback(hObject, eventdata, handles)
global im_file

[file1,path1]=uigetfile({'*.tif;*.tiff','TIF Image Files'},'Select Signal Stack');
if file1==0
    return
end
cd(path1)

im_file=fullfile(path1,file1);

set(handles.text_status,'String',['loaded: ' file1])


% --- Executes on button press in pushbutton_analyze.
function pushbutton_analyze_Callback(hObject, eventdata, handles)
global mat_file
global im_file

global dl
global res
global offset

global rel_sz
global rel_sigma

span=str2num(get(handles.edit_span,'String'));
degree=str2num(get(handles.edit_degree,'String'));

if or(isempty(mat_file),isempty(im_file))
    warndlg('Load a *.mat and image stack file to begin analysis.','Load files first.')
    return
end

load(mat_file)
Nframe=length(frame);
Ncell=0;
for i=1:Nframe
    temp1=max([frame(i).object.cellID]);
    if temp1>Ncell
        Ncell=temp1;
    end
end
set(handles.text_status,'String',['Found ' num2str(Nframe) ' frames and ' num2str(Ncell) 'cells.'])

%get shifts
dx=str2num(get(handles.edit_xshift,'String'));
dy=str2num(get(handles.edit_yshift,'String'));

%check cell ID questions
all_cells=get(handles.checkbox_allcells,'Value');
ignore1=get(handles.checkbox_ignoreproximal,'Value');
cell_num=str2num(get(handles.edit_cellnumber,'String'));

contq=get(handles.checkbox_contour,'Value');
intq=get(handles.checkbox_interior,'Value');
pillq=get(handles.checkbox_pill,'Value');

%check for pill mesh data
if pillq
    if ~isfield(frame(1).object(1),'pill_mesh')
        if get(handles.checkbox_contour,'Value')
            warndlg('This mat-file has no pill mesh data. Processing will continue on contours only.')
        else
            warndlg('This mat-file has no pill mesh data. Processing has stoppped.')
            return
        end
        set(handles.checkbox_pill,'Value',0)
        set(handles.checkbox_normalize,'Enable','off')
    end
end

%check for cell table
if get(handles.checkbox_singlecell,'Value')
    if ~exist('cells','var')
        set(handles.text_status,'String',['Creating Cell ID Table.'])
        %check data types
        if v_indpt==0
            chk1=1;
        else
            chk1=0;
            warning('The cells in this file may have been processed as separate objects.')
        end
        
        if isfield(frame(1).object,'cellID')
            chk2=1;
        else
            chk2=0;
            error('Could not find the necessary data to construct the cell table.')
        end
        
        %construct table
        for i=1:Ncell
            set(handles.text_status,'String',['Processing cell: ' num2str(i)])
            q=0;
            for j=1:length(frame)
                for k=1:frame(j).num_objs
                    ind1=frame(j).object(k).cellID;
                    
                    if ind1==i
                        q=q+1;
                        cells(i).frame(q)=j;
                        cells(i).object(q)=k;
                    end
                end
            end
        end
        
        %save table
        save(mat_file,'cells','-append')
        
        set(handles.text_status,'String',['Cell tabled created.'])
    end
    
    if length(cells)<cell_num
        warndlg(['You requested processing cell ' num2str(cell_num) ', but the highest found cell number was ' num2str(length(cells)) '.'],'Cell Number Error')
        return
    end
end

%get channel name
prompt={'What is the name of this channel? (e.g. BF, GFP)'};
name='Input channel name';
numlines=1;
defaultanswer={''};

answer1=inputdlg(prompt,name,numlines,defaultanswer);

if isempty(answer1)
    return
elseif isempty(answer1{1})
    return
end


%****************************************
%****** BACKGROUND **********************
%****************************************

%calculate background
if get(handles.checkbox_background,'Value')
    %check current background types in mat-file
    if isfield(frame,'background')
        %get previous channel names
        nameset={};
        for i=1:length(frame(1).background)
            nameset{i}=frame(1).background(i).name;
        end
        
        if any(strcmp(answer1{1},nameset))
            %ask to overwrite
            quest1 = questdlg(['A background channel named ' answer1{1} ' already exists; overwrite this data?'], ...
                'Overwrite?','Yes', 'No', 'Yes');
            
            if strcmp(quest1,'No')
                return
            else
                n=find(strcmp(answer1{1},nameset));
            end
        else
            n=length(nameset)+1;
        end
    else
        n=1;
    end
    
    plotq=get(handles.checkbox_plotq,'Value');
    
    %create master background image
    T(:,:,1)=mat2gray(imread(im_file,1));
    T(:,:,2)=mat2gray(imread(im_file,round(Nframe/2)));
    T(:,:,3)=mat2gray(imread(im_file,Nframe));
    
    f1=figure;
    imagesc(mat2gray(T))
    axis equal tight
    title(['Select Representative Background Region'])
    xlabel('X')
    ylabel('Y')
    set(gcf, 'Position', get(0,'Screensize'));
    
    %clear rect
    rect1=zeros(1,4);
    
    while (rect1(3)*rect1(4))<=25
        rect1=getrect;
        
        %calculate box coords
        xmin=round(rect1(1));
        xmax=round(rect1(1)+rect1(3));
        ymin=round(rect1(2));
        ymax=round(rect1(2)+rect1(4));
    end
    
    %process multiple channels
    back_mean=zeros(Nframe,1);
    back_median=zeros(Nframe,1);
    back_std=zeros(Nframe,1);
    
    for i=1:Nframe
        Im1=double(imread(im_file,i));
        
        %apply relnoise filter
        if get(handles.checkbox_noise,'Value')
            Im1=relnoise(Im1,rel_sz,rel_sigma);
        end
        
        %region stats
        region1=double(Im1(ymin:ymax,xmin:xmax));
        back_mean(i)=mean(region1(:));
        back_median(i)=median(region1(:));
        back_std(i)=std(region1(:));
        
        %store values in 'frame' structure
        frame(i).background(n).name=answer1{1};
        frame(i).background(n).mean=back_mean(i);
        frame(i).background(n).median=back_median(i);
        frame(i).background(n).std=back_std(i);
        
        if plotq
            %update figure
            hold off
            imagesc(Im1)
            hold on
            colormap(gray)
            plot([xmin,xmax],[ymin,ymin],'w--')
            plot([xmin,xmax],[ymax,ymax],'w--')
            plot([xmin,xmin],[ymin,ymax],'w--')
            plot([xmax,xmax],[ymin,ymax],'w--')
            title(['Frame ' num2str(i) ' of ' num2str(Nframe)])
            axis equal tight
            drawnow
        else
            %clc
            title(['Frame ' num2str(i) ' of ' num2str(Nframe)])
            drawnow
        end
    end
    close(f1)
    
    set(handles.text_status,'String',['Background data will be added to mat-file.'])
end

%check if kymograph data is already present
if isfield(frame(1).object(1),'channel')
    nameset={};
    for i=1:length(frame(1).object(1).channel)
        nameset{i}=frame(1).object(1).channel(i).name;
    end
    
    if any(strcmp(answer1{1},nameset))
        %ask to overwrite
        quest1 = questdlg(['A signal channel named ' answer1{1} ' already exists; overwrite this data?'], ...
            'Overwrite?','Yes', 'No', 'Yes');
        
        if strcmp(quest1,'No')
            return
        else
            m=find(strcmp(answer1{1},nameset));
        end
    else
        m=length(nameset)+1;
    end
else
    m=1;
end

tic
if get(handles.checkbox_plotq,'Value')
    f1=figure;
end

%frame loop
for i=1:Nframe
    %get cell list
    if all_cells
        if ignore1
            %find proximal cells
            prox1=unique([frame(i).object.proxID]);
        else
            prox1=[];
        end
        
        cell_list=setdiff(1:length(frame(i).object),prox1);
    else
        ind1=find(cells(cell_num).frame==i,1,'first');
        if ~isempty(ind1)
            cell_list=cells(cell_num).object(ind1);
        else
            cell_list=[];
        end
    end
    
    if isempty(cell_list)
        disp(['Frame ' num2str(i) ' had no valid cells to process.'])
        continue
    end
    
    %apply relnoise filter
    im1=imread(im_file,i);
    if get(handles.checkbox_noise,'Value')
        im1=relnoise(im1,rel_sz,rel_sigma);
    end
    
    if get(handles.checkbox_plotq,'Value')
        figure(f1);
        subplot(1,2,1)
        hold off
        imagesc(im1)
        hold on
        xlabel('X')
        ylabel('Y')
        title(['Frame ' num2str(i)])
        box on
        axis equal tight
        colormap gray
        
        %figure(f1);
        %subplot(1,2,2)
        %cla
    end
    
    for j=cell_list
        X=dx+frame(i).object(j).Xcont;
        Y=dy+frame(i).object(j).Ycont;
        
        %set channel name
        frame(i).object(j).channel(m).name=answer1{1};
        
        %set background type
        if and(get(handles.checkbox_background,'Value'),get(handles.checkbox_subtract,'Value'))
            frame(i).object(j).channel(m).back_type='subtracted';
        else
            frame(i).object(j).channel(m).back_type='none';
        end
        
        %*********************************
        %**** CONTOUR INTENSITY **********
        %*********************************
        if contq
            prof=get_contour_profile(X,Y,im1,dl,res,offset);
            
            %polynomial smooth
            if get(handles.checkbox_smooth,'Value')
                poly_temp=smooth((1:(3*length(X)))',[prof;prof;prof],span,'sgolay',degree);
                prof=poly_temp(length(X)+1:(2*length(X)));
            end
            
            %subtract background
            if and(get(handles.checkbox_background,'Value'),get(handles.checkbox_subtract,'Value'))
                prof=prof-back_mean(i);
            end
            
            frame(i).object(j).channel(m).contour_int=prof;
        end
        
        
        %**********************************
        %**** INTERIOR INTENSITY **********
        %**********************************
        if intq
            %perimeter pixels
            x1=round(frame(i).object(j).Xperim)+dx;
            y1=round(frame(i).object(j).Yperim)+dy;
            
            %reduced image size
            sx_tiny=max(x1)-min(x1)+1;
            sy_tiny=max(y1)-min(y1)+1;
            
            %reduced image coord shift
            dx_tiny=x1-min(x1)+1;
            dy_tiny=y1-min(y1)+1;
            
            %create sub image from original
            im1_tiny=im1(min(y1):min(y1)+sy_tiny-1,min(x1):min(x1)+sx_tiny-1);
            
            tiny_perim=zeros(sy_tiny,sx_tiny);
            for k=1:length(x1)
               tiny_perim(dy_tiny(k),dx_tiny(k))=1;
            end
            
            %create label image
            tiny_mask=imfill(tiny_perim,'holes');
            area1=sum(tiny_mask(:));
            
            %get data
            S2=regionprops(tiny_mask,im1_tiny,'PixelValues','MeanIntensity');
            mean1=[S2.MeanIntensity];
            
            %subtract background
            if and(get(handles.checkbox_background,'Value'),get(handles.checkbox_subtract,'Value'))
                mean1=mean1-back_mean(i);
            end
            
            frame(i).object(j).channel(m).mean_int=mean1; %This is the internal mean intensity
            frame(i).object(j).channel(m).mean_int_area=area1;
        end
        
        %*******************************************
        %**** PILL MESH INTERNAL INTENSITY *********
        %*******************************************
        if pillq
            %create polygons from mesh
            n_areas=size(frame(i).object(j).pill_mesh,1)-1;
            
            xpoly=zeros(n_areas-1,4);
            ypoly=zeros(n_areas-1,4);
            for k=1:n_areas
                %xpoly(k,:)=[frame(i).object(j).pill_mesh(k,1),frame(i).object(j).pill_mesh(k,3),frame(i).object(j).pill_mesh(k+1,3),frame(i).object(j).pill_mesh(k+1,1)];
                %ypoly(k,:)=[frame(i).object(j).pill_mesh(k,2),frame(i).object(j).pill_mesh(k,4),frame(i).object(j).pill_mesh(k+1,4),frame(i).object(j).pill_mesh(k+1,2)];
                
                xpoly(k,:)=[frame(i).object(j).pill_mesh(k+1,1),frame(i).object(j).pill_mesh(k+1,3),frame(i).object(j).pill_mesh(k,3),frame(i).object(j).pill_mesh(k,1)];
                ypoly(k,:)=[frame(i).object(j).pill_mesh(k+1,2),frame(i).object(j).pill_mesh(k+1,4),frame(i).object(j).pill_mesh(k,4),frame(i).object(j).pill_mesh(k,2)];
            end
            
            pill_int=zeros(n_areas,1);
            for k=1:n_areas
                pill_int(k)=get_pillmesh_profile(xpoly(k,:),ypoly(k,:),im1);
                set(handles.text_status,'String',['Computing vertex ' num2str(k) ', object ' num2str(j) ', frame ' num2str(i) '.'])
                drawnow
            end
            
            if get(handles.checkbox_normalize,'Value')
                %calculate mesh quadrangle areas
                pill_area=polyarea(xpoly',ypoly')';
                
                %normalize integrated intensity
                pill_int=pill_int./pill_area; 
                
                %subtract background
                if and(get(handles.checkbox_background,'Value'),get(handles.checkbox_subtract,'Value'))
                    pill_int=pill_int-back_mean(i);
                end
                
                %set normalization type
                %Storing variables problem was fixed. Storing pill_area
                %only when normalizing
                % Andres Florez - 05/17/21 
                frame(i).object(j).channel(m).internal_norm='normalized';
                frame(i).object(j).channel(m).internal_int_norm=pill_int;
                frame(i).object(j).channel(m).pill_mesh_area=pill_area;    
            else
                frame(i).object(j).channel(m).internal_norm='unnormalized';
            end
            %if no normalization is done it won't save pill area since it
            %is only calculated in the previous step - Andres Florez
            %05/17/21
            frame(i).object(j).channel(m).internal_int=pill_int;
            
            %{
            %test figure
            figure;
            hold on
            plot(Xbox,Ybox,'b-')
            plot(Xpoly',Ypoly','r-')
            %}
        end
        %*******************************************
        
        %initialize figure
        if get(handles.checkbox_plotq,'Value')
            %plot contour
            figure(f1);
            subplot(1,2,1)
            plot(X,Y,'r')
            text(mean(X),mean(Y),num2str(frame(i).object(j).cellID),'color',[0,0.9,1])
            
            figure(f1);
            subplot(1,2,2)
            if and(contq,~pillq)
                plot(prof,'b','Linewidth',2)
                ylabel('Intensity')
                xlabel('Contour Point')
                xlim([1 length(X)])
                title(['Contour intensity profile for object ' num2str(j)])
                
            elseif and(~contq,pillq)
                plot(pill_int,'r','Linewidth',2)
                if get(handles.checkbox_normalize,'Value')
                    ylabel('Intensity (area normalized)')
                else
                    ylabel('Intensity')
                end
                xlabel('Centerline Point')
                xlim([1 length(pill_int)])
                title(['Centerline intensity profile for object ' num2str(j)])
            elseif and(contq,pillq) % Added this condition to fix the plotting - Andres Florez 05/17/21
                hold off
                plot(prof,'b','Linewidth',2)
                hold on
                plot(pill_int,'r','Linewidth',2)
                
                if get(handles.checkbox_normalize,'Value')
                    ylabel('Intensity (area normalized)')
                else
                    ylabel('Intensity')
                end
                xlabel('Centerline Point')
                xlim([1 length(X)])
                title(['Centerline (r) and contour (b) intensity profiles for object ' num2str(j)])
            end
            box on
            drawnow
        else
            set(handles.text_status,'String',['Processing object ' num2str(j) ' on frame ' num2str(i) '.'])
        end
    end
end

%perform single cell contour alignment
if all([contq,~all_cells,get(handles.checkbox_align,'Value')])
    
end

delete(mat_file)
try
    save(mat_file,'f_*','v_*','frame','cells')
catch
    save(mat_file,'f_*','v_*','frame')
end

set(handles.text_status,'String',['Finished processing in ' num2str(round(100*toc/60)/100) ' minutes.'])









function edit_cellnumber_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_cellnumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_xshift_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_xshift_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_yshift_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_yshift_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_background.
function checkbox_background_Callback(hObject, eventdata, handles)
if get(handles.checkbox_background,'Value')
    set(handles.checkbox_subtract,'Enable','on')
else
    set(handles.checkbox_subtract,'Enable','off')
end

% --- Executes on button press in checkbox_ignoreproximal.
function checkbox_ignoreproximal_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_contour.
function checkbox_contour_Callback(hObject, eventdata, handles)
if get(handles.checkbox_contour,'Value')
    set(handles.checkbox_align,'Enable','on')
    set(handles.checkbox_smooth,'Enable','on')
    set(handles.edit_span,'Enable','on')
    set(handles.edit_degree,'Enable','on')
    set(handles.text5,'Enable','on')
    set(handles.text6,'Enable','on')
else
    set(handles.checkbox_align,'Enable','off')
    set(handles.checkbox_smooth,'Enable','off')
    set(handles.edit_span,'Enable','off')
    set(handles.edit_degree,'Enable','off')
    set(handles.text5,'Enable','off')
    set(handles.text6,'Enable','off')
end

if all([~get(handles.checkbox_pill,'Value'),~get(handles.checkbox_contour,'Value'),~get(handles.checkbox_interior,'Value')])
    set(handles.pushbutton_analyze,'Enable','off')
    set(handles.checkbox_plotq,'Enable','off')
    set(handles.checkbox_noise,'Enable','off')
else
    set(handles.pushbutton_analyze,'Enable','on')
    set(handles.checkbox_plotq,'Enable','on')
    set(handles.checkbox_noise,'Enable','on')
end

% --- Executes on button press in checkbox_interior.
function checkbox_interior_Callback(hObject, eventdata, handles)
if all([~get(handles.checkbox_pill,'Value'),~get(handles.checkbox_contour,'Value'),~get(handles.checkbox_interior,'Value')])
    set(handles.pushbutton_analyze,'Enable','off')
    set(handles.checkbox_plotq,'Enable','off')
    set(handles.checkbox_noise,'Enable','off')
else
    set(handles.pushbutton_analyze,'Enable','on')
    set(handles.checkbox_plotq,'Enable','on')
    set(handles.checkbox_noise,'Enable','on')
end

% --- Executes on button press in checkbox_pill.
function checkbox_pill_Callback(hObject, eventdata, handles)

if ~license('test','map_toolbox')
    warndlg('Your version of Matlab(c) seems to be missing the Mapping Toolbox required for this functionality.','Missing Toolbox');
    set(handles.checkbox_pill,'Value',0)
    return
end

if get(handles.checkbox_pill,'Value')
    set(handles.checkbox_normalize,'Enable','on')
    set(handles.text_status,'String','This analysis can take a long time.')
else
    set(handles.checkbox_normalize,'Enable','off')
    set(handles.text_status,'String',' ')
end

if all([~get(handles.checkbox_pill,'Value'),~get(handles.checkbox_contour,'Value'),~get(handles.checkbox_interior,'Value')])
    set(handles.pushbutton_analyze,'Enable','off')
    set(handles.checkbox_plotq,'Enable','off')
    set(handles.checkbox_noise,'Enable','off')
else
    set(handles.pushbutton_analyze,'Enable','on')
    set(handles.checkbox_plotq,'Enable','on')
    set(handles.checkbox_noise,'Enable','on')
end

% --- Executes on button press in checkbox_subtract.
function checkbox_subtract_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_align.
function checkbox_align_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_smooth.
function checkbox_smooth_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_normalize.
function checkbox_normalize_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_noise.
function checkbox_noise_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_allcells.
function checkbox_allcells_Callback(hObject, eventdata, handles)
set(handles.checkbox_allcells,'Value',1)
set(handles.checkbox_ignoreproximal,'Enable','on')
set(handles.checkbox_singlecell,'Value',0)
set(handles.edit_cellnumber,'Enable','off')

% --- Executes on button press in checkbox_singlecell.
function checkbox_singlecell_Callback(hObject, eventdata, handles)
set(handles.checkbox_singlecell,'Value',1)
set(handles.checkbox_ignoreproximal,'Enable','off')
set(handles.checkbox_allcells,'Value',0)
set(handles.edit_cellnumber,'Enable','on')

% --- Executes on button press in checkbox_plotq.
function checkbox_plotq_Callback(hObject, eventdata, handles)



function edit_span_Callback(hObject, eventdata, handles)
temp1=str2num(get(handles.edit_span,'String'));
temp2=str2num(get(handles.edit_degree,'String'));

if temp1<4
    temp1=4;
end

if temp1<=(temp2+2)
    temp1=temp2+3;
end
set(handles.edit_span,'String',num2str(temp1))

% --- Executes during object creation, after setting all properties.
function edit_span_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_degree_Callback(hObject, eventdata, handles)
temp1=str2num(get(handles.edit_span,'String'));
temp2=str2num(get(handles.edit_degree,'String'));

if temp2<3
    temp1=3;
end

if temp1<=(temp2+2)
    temp2=temp1-3;
end
set(handles.edit_degree,'String',num2str(temp2))

% --- Executes during object creation, after setting all properties.
function edit_degree_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function options_1_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function contour_parameters_Callback(hObject, eventdata, handles)
contour_parameters

% --------------------------------------------------------------------
function noise_filt_parameters_Callback(hObject, eventdata, handles)
noise_filt_parameters

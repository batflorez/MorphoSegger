function varargout = quality_check(varargin)
% QUALITY_CHECK MATLAB code for quality_check.fig
%      QUALITY_CHECK, by itself, creates a new QUALITY_CHECK or raises the existing
%      singleton*.
%
%      H = QUALITY_CHECK returns the handle to a new QUALITY_CHECK or the handle to
%      the existing singleton*.
%
%      QUALITY_CHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUALITY_CHECK.M with the given input arguments.
%
%      QUALITY_CHECK('Property','Value',...) creates a new QUALITY_CHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quality_check_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quality_check_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quality_check

% Last Modified by GUIDE v2.5 17-Aug-2015 00:51:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @quality_check_OpeningFcn, ...
    'gui_OutputFcn',  @quality_check_OutputFcn, ...
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


% --- Executes just before quality_check is made visible.
function quality_check_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for quality_check
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

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = quality_check_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in pushbutton_qualitycheck.
function pushbutton_qualitycheck_Callback(hObject, eventdata, handles)
global mat_file
global im_file

if or(isempty(mat_file),isempty(im_file))
    warndlg('Load a *.mat and image stack file to begin quality check.','Load files first.')
    return
end

load(mat_file)

Nframe = length(frame);
Ncell=max([frame(end).object.cellID]);
set(handles.text_status,'String',['Found ' num2str(Nframe) ' frames and ' num2str(Ncell) 'cells.'])

%get shifts
dx=str2num(get(handles.edit_xshift,'String'));
dy=str2num(get(handles.edit_yshift,'String'));

%check cell ID questions
all_cells=get(handles.checkbox_allcells,'Value');
ignore1=get(handles.checkbox_ignoreproximal,'Value');
cell_num=str2num(get(handles.edit_cellnumber,'String'));

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


tic
f1=figure;

sz1=size(imread(im_file,1));
break_opt=0;

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
    
    im1=imread(im_file,i);
    
    for j=cell_list
        X=dx+frame(i).object(j).Xcont;
        Y=dy+frame(i).object(j).Ycont;
        
        %get bounds
        buffer1=10; %pxls
        
        minX=round(min(X(:))-buffer1);
        maxX=round(max(X(:))+buffer1);
        minY=round(min(Y(:))-buffer1);
        maxY=round(max(Y(:))+buffer1);
        
        minX(minX<1)=1;
        maxX(maxX>sz1(2))=sz1(2);
        minY(minY<1)=1;
        maxY(maxY>sz1(1))=sz1(1);
        
        %clip image
        im1_sub=im1(minY:maxY,minX:maxX);
        
        %reposition
        X=X-(minX-1);
        Y=Y-(minY-1);
        
        %plot contour
        figure(f1);
        imagesc(im1_sub)
        hold on
        plot(X,Y,'r','linewidth',2)
        hold off
        xlabel('X')
        ylabel('Y')
        title(['Frame ' num2str(i) ', obejct ' num2str(j)])
        box on
        axis equal tight
        colormap gray
        drawnow
        
        cont1=1;
        while cont1
            %get user input
            %select point in current region of interest
            button1=getkey;
            
            %if user clicks on non-region, give option list
            if ~isempty(button1)
                if sum(button1==[31,29,8])>0
                    %down arrow = goes forward a frame and rejects contour
                    %right arrow = goes forward a frame and approves contour
                    %delete/backspace = halt quality check on all images
                    
                    cont1=0;
                    
                    %right arrow
                    if button1==29
                        frame(i).object(j).quality=1;
                    end
                    
                    %down arrow
                    if button1==31
                        frame(i).object(j).quality=0;
                    end
                    
                    %delete
                    if button1==8
                        break_opt=1;
                        break
                    end
                else
                    title('Invalid key stroke.')
                    pause(1)
                    title(['Frame ' num2str(i) ', obejct ' num2str(j)])
                end
            else
                title('Invalid key stroke.')
                pause(1)
                title(['Frame ' num2str(i) ', obejct ' num2str(j)])
            end
        end
        
        
        if break_opt
            break
        end
    end
    
    if break_opt
        break
    end
end
close(f1)

%modify frame data structuredelete(mat_file)
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


% --- Executes on button press in checkbox_ignoreproximal.
function checkbox_ignoreproximal_Callback(hObject, eventdata, handles)

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

function varargout = stack_align(varargin)
% STACK_ALIGN MATLAB code for stack_align.fig
%      STACK_ALIGN, by itself, creates a new STACK_ALIGN or raises the existing
%      singleton*.
%
%      H = STACK_ALIGN returns the handle to a new STACK_ALIGN or the handle to
%      the existing singleton*.
%
%      STACK_ALIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_ALIGN.M with the given input arguments.
%
%      STACK_ALIGN('Property','Value',...) creates a new STACK_ALIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_align_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_align_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_align

% Last Modified by GUIDE v2.5 02-Sep-2015 13:07:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_align_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_align_OutputFcn, ...
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


% --- Executes just before stack_align is made visible.
function stack_align_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_align (see VARARGIN)

% Choose default command line output for stack_align
handles.output = hObject;

%input file info
set(handles.edit_namemotif,'String','Pos*.tif') %name motif
set(handles.edit_channelmotif,'String','_YFP, _mCherry') %channel motif

%align stacks
set(handles.checkbox_master1,'Value',1) %set the master stack for alignment
set(handles.checkbox_master2,'Value',0) %set the master stack for alignment
set(handles.checkbox_master3,'Value',0) %set the master stack for alignment
set(handles.checkbox_master4,'Value',0) %set the master stack for alignment

set(handles.checkbox_slave1,'Value',0) %set the slave stack for alignment
set(handles.checkbox_slave2,'Value',0) %set the slave stack for alignment
set(handles.checkbox_slave3,'Value',0) %set the slave stack for alignment
set(handles.checkbox_slave4,'Value',0) %set the slave stack for alignment

set(handles.edit_alignshift,'String','10') %set maximum XY shift for alignment

edit_channelmotif_Callback(hObject, eventdata, handles)
checkbox_master1_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = stack_align_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_namemotif_Callback(hObject, eventdata, handles)
if isempty(get(handles.edit_namemotif,'String'))
    set(handles.edit_namemotif,'String','*.tif')
end

%count the number of files matching the motif
dir0=dir(fullfile(cd,get(handles.edit_namemotif,'String')));
if isempty(dir0)
    nfiles=0;
else
    nfiles=sum([dir0.isdir]==0);
end
set(handles.text_filesfound,'String',[num2str(nfiles) ' files match the image motif.'])

% --- Executes during object creation, after setting all properties.
function edit_namemotif_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_channelmotif_Callback(hObject, eventdata, handles)
%parse inputs
channelname=parse_input(get(handles.edit_channelmotif,'String'));

if length(channelname)>4
    set(handles.edit_channelmotif,'String',[channelname{1} ', ' channelname{2} ', ' channelname{3} ', ' channelname{4}])
    errordlg('This feature currently supports up to four input channels.', 'Too many channels');
    return
elseif isempty(channelname)
    set(handles.edit_channelmotif,'String','enter a channel')
    errordlg('You must specify at least one, and no more than four channels.');
    return
else
    for i=1:4
        %set text
        if i<=length(channelname)
            temp_name=channelname{i};
            eval(['set(handles.checkbox_master' num2str(i) ',''String'',temp_name(temp_name~=''_''))'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''String'',temp_name(temp_name~=''_''))'])
        else
            eval(['set(handles.checkbox_master' num2str(i) ',''String'',[''Ch'' num2str(i) ])'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''String'',[''Ch'' num2str(i) ])'])
        end
        
        %set active GUI regions
        if i<=length(channelname)
            eval(['set(handles.checkbox_master' num2str(i) ',''Enable'',''On'')'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''Enable'',''On'')'])
        
            eval(['set(handles.checkbox_master' num2str(i) ',''Visible'',''On'')'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''Visible'',''On'')'])
        else
            eval(['set(handles.checkbox_master' num2str(i) ',''Enable'',''Off'')'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''Enable'',''Off'')'])
            
            eval(['set(handles.checkbox_master' num2str(i) ',''Visible'',''Off'')'])
            eval(['set(handles.checkbox_slave' num2str(i) ',''Visible'',''Off'')'])
        end
    end
end

switch 1
    case get(handles.checkbox_master1,'Value')
        set(handles.checkbox_slave1,'Value',0)
        set(handles.checkbox_slave1,'enable','off')
    case get(handles.checkbox_master2,'Value')
        set(handles.checkbox_slave2,'Value',0)
        set(handles.checkbox_slave2,'enable','off')
    case get(handles.checkbox_master3,'Value')
        set(handles.checkbox_slave3,'Value',0)
        set(handles.checkbox_slave3,'enable','off')
    case get(handles.checkbox_master4,'Value')
        set(handles.checkbox_slave4,'Value',0)
        set(handles.checkbox_slave4,'enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit_channelmotif_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_alignshift_Callback(hObject, eventdata, handles)
temp1=str2num(get(handles.edit_alignshift,'String'));
if temp1<1
    set(handles.edit_alignshift,'String','1')
elseif temp1>49
    warndlg('Keep in mind that the time required to align images grows quadratically with this number.')
end
% --- Executes during object creation, after setting all properties.
function edit_alignshift_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_master1.
function checkbox_master1_Callback(hObject, eventdata, handles)
set(handles.checkbox_master1,'Value',1)
set(handles.checkbox_master2,'Value',0)
set(handles.checkbox_master3,'Value',0)
set(handles.checkbox_master4,'Value',0)

set(handles.checkbox_slave1,'Value',0)
set(handles.checkbox_slave1,'enable','off')
set(handles.checkbox_slave2,'enable','on')
set(handles.checkbox_slave3,'enable','on')
set(handles.checkbox_slave4,'enable','on')

% --- Executes on button press in checkbox_master2.
function checkbox_master2_Callback(hObject, eventdata, handles)
set(handles.checkbox_master1,'Value',0)
set(handles.checkbox_master2,'Value',1)
set(handles.checkbox_master3,'Value',0)
set(handles.checkbox_master4,'Value',0)

set(handles.checkbox_slave2,'Value',0)
set(handles.checkbox_slave1,'enable','on')
set(handles.checkbox_slave2,'enable','off')
set(handles.checkbox_slave3,'enable','on')
set(handles.checkbox_slave4,'enable','on')

% --- Executes on button press in checkbox_master3.
function checkbox_master3_Callback(hObject, eventdata, handles)
set(handles.checkbox_master1,'Value',0)
set(handles.checkbox_master2,'Value',0)
set(handles.checkbox_master3,'Value',1)
set(handles.checkbox_master4,'Value',0)

set(handles.checkbox_slave3,'Value',0)
set(handles.checkbox_slave1,'enable','on')
set(handles.checkbox_slave2,'enable','on')
set(handles.checkbox_slave3,'enable','off')
set(handles.checkbox_slave4,'enable','on')

% --- Executes on button press in checkbox_master4.
function checkbox_master4_Callback(hObject, eventdata, handles)
set(handles.checkbox_master1,'Value',0)
set(handles.checkbox_master2,'Value',0)
set(handles.checkbox_master3,'Value',0)
set(handles.checkbox_master4,'Value',1)

set(handles.checkbox_slave4,'Value',0)
set(handles.checkbox_slave1,'enable','on')
set(handles.checkbox_slave2,'enable','on')
set(handles.checkbox_slave3,'enable','on')
set(handles.checkbox_slave4,'enable','off')

% --- Executes on button press in checkbox_slave1.
function checkbox_slave1_Callback(hObject, eventdata, handles)
if get(handles.checkbox_master1,'Value')
    set(handles.checkbox_slave1,'Value',0)
end

% --- Executes on button press in checkbox_slave2.
function checkbox_slave2_Callback(hObject, eventdata, handles)
if get(handles.checkbox_master2,'Value')
    set(handles.checkbox_slave2,'Value',0)
end

% --- Executes on button press in checkbox_slave3.
function checkbox_slave3_Callback(hObject, eventdata, handles)
if get(handles.checkbox_master3,'Value')
    set(handles.checkbox_slave3,'Value',0)
end

% --- Executes on button press in checkbox_slave4.
function checkbox_slave4_Callback(hObject, eventdata, handles)
if get(handles.checkbox_master4,'Value')
    set(handles.checkbox_slave4,'Value',0)
end

% --- Executes on button press in checkbox_adaptive.
function checkbox_adaptive_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_alignimages.
function pushbutton_alignimages_Callback(hObject, eventdata, handles)
%Tristan Ursell
%Align Images by Contour
%Feb 2013
%
% align_images
%
%This script uses a cell's contour to generate a contour-centroid-centered
%set of images from the input images, i.e. to get rid of drift (but not
%rotation, scaling, or skew).
%
%The input parmateter 'dx_max' sets the size of the local translation
%window about with the cross-correlation will be calculated.  Specifically,
%the correlation matrix will of size (2*dx_max+1) x (2*dx_max+1).
%
%The input parameter 'filepath' is a string the defining the full file name
%and path.  If not 'filepath' is specificed, a GUI window will pop up for
%file selection.

%maximum shift in either x or y in between frames
dx_max=str2num(get(handles.edit_alignshift,'String'));
dx_max0=dx_max;

%get channels
channelname=parse_input(get(handles.edit_channelmotif,'String'));

%get basal name
basename=strtrim(get(handles.edit_namemotif,'String'));

[file0,path0,~]=uigetfile({'*.mat;*.tif;*.tiff','mat & TIF files'},'Directory selection -- pick any file',basename);
if file0==0
    return
end
cd(path0)

%creat basal file list without directories listed
temp1=dir(basename);
dir_base=temp1(~[temp1.isdir]);

if isempty(dir_base)
    errordlg(['No files in this directory matched the search terms.'],'No files found');
    set(handles.text_filesfound,'String',['0 files match the image motif.'])
    return
end

%get master and slave strings
for i=1:length(channelname)
    master_vec(i)=eval(['get(handles.checkbox_master' num2str(i) ',''Value'')']);
    slave_vec(i)=eval(['get(handles.checkbox_slave' num2str(i) ',''Value'')']);
end

%label master
types_vec=zeros(length(dir_base),1);
for i=1:length(dir_base)
    %find master
    if ~isempty(strfind(dir_base(i).name,channelname{find(master_vec)}))
        types_vec(i)=-1;
    end
    
    %find slaves
    for j=find(slave_vec)
        if ~isempty(strfind(dir_base(i).name,channelname{j}))
            types_vec(i)=j;
        end
    end
end

%check for too many masters
if sum(types_vec==-1)>1
    errordlg('Found too many master files; try making ''Image motif'' more specific or did you forget to combine images into a stack?')
    return
end

%get master file 
master_file=dir_base(types_vec==-1).name;

%get image info
info1=imfinfo(master_file);
sz_x=info1.Width;
sz_y=info1.Height;
N=length(info1);

%select calculation region
clear T

temp1=imread(master_file,1);
temp2=imread(master_file,round(N/2));
temp3=imread(master_file,N);

T(:,:,1)=mat2gray(temp1);
T(:,:,2)=mat2gray(temp2);
T(:,:,3)=mat2gray(temp3);

%select alignment window from first frame

f1=figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(T)
axis equal tight
title('Select Alignment Region')
rect1=round(getrect);

%initialize shift vectors
X_shift=zeros(N,1);
Y_shift=zeros(N,1);

cptQ = questdlg('Choose a center point?','Center Point','Yes', 'No', 'No');

if strcmp(cptQ,'Yes')
    title('Select Centering Point for Output Images (''backspace'' to delete current point, ''enter'' to finish)')
    [pt_x,pt_y]=getpts;
    
    X_shift(1)=round(sz_x/2-pt_x);
    Y_shift(1)=round(sz_y/2-pt_y);
end

%{
%align to region center
X_shift(1)=round(sz_x/2-(rect1(1)+rect1(3)/2));
Y_shift(1)=round(sz_y/2-(rect1(2)+rect1(4)/2));
%}

close(f1);
drawnow

%determine shift vector for each frame
for i=2:N
    %get images
    curr_im1=double(imread(master_file,i-1));
    curr_im1=circshift(curr_im1,[Y_shift(i-1),X_shift(i-1)]);
    curr_im2=double(imread(master_file,i));
    
    %get sub-image and normalize
    sub_im1_raw=curr_im1(rect1(2):rect1(2)+rect1(4),rect1(1):rect1(1)+rect1(3));
    sub_im1_mean=sub_im1_raw-mean(sub_im1_raw(:));
    sub_im1=sub_im1_mean/sqrt(sum(sub_im1_mean(:).^2));
    
    %reset cross corr matrix
    cc_mat1=zeros(2*dx_max0+1,2*dx_max0+1);
    
    %calculate the cross-correlation (or similar) matrix
    dy=0;
    ind_vecy=Y_shift(i-1)+(-dx_max0:dx_max0);
    ind_vecx=X_shift(i-1)+(-dx_max0:dx_max0);
    for q1=ind_vecy
        %perform row shift
        dy=dy+1;
        curr_im2_temp=circshift(curr_im2,[q1,0]);
        curr_im2_temp_crop=curr_im2_temp(rect1(2):rect1(2)+rect1(4),:);
        
        dx=0;
        for q2=ind_vecx
            %perform column shift
            dx=dx+1;
            curr_im2_temp2=circshift(curr_im2_temp_crop,[0,q2]);
            
            %clip out correct region and normalize
            sub_im2_raw=curr_im2_temp2(:,rect1(1):rect1(1)+rect1(3));
            sub_im2_mean=sub_im2_raw-mean(sub_im2_raw(:));
            sub_im2=sub_im2_mean/sqrt(sum(sub_im2_mean(:).^2));
            
            %calculate correlation coefficient
            cc_mat1(dy,dx)=sum(sub_im1(:).*sub_im2(:));
        end
            set(handles.text_progress,'String',['computing shifts: ' num2str(round(100*i/N)) '%'])
            drawnow
    end
    
    %find shift values
    [y0,x0]=find(cc_mat1==max(cc_mat1(:)),1,'first');
    
    %check if on edge
    if or(or(y0==1,y0==2*dx_max0+1),or(x0==1,x0==2*dx_max0+1))
        disp(['Warning: shift limit reached on frame ' num2str(i) '.'])
        bnd_warn=1;
    else
        bnd_warn=0;
    end
    
    %readjust if boundary of region reached
    dx_max=dx_max0;
    while and(bnd_warn,get(handles.checkbox_adaptive,'Value'))
        dx_max=dx_max+2;
        %reset cross corr matrix
        cc_mat1=zeros(2*dx_max+1,2*dx_max+1);
        
        %calculate the cross-correlation (or similar) matrix
        dy=0;
        ind_vecy=Y_shift(i-1)+(-dx_max:dx_max);
        ind_vecx=X_shift(i-1)+(-dx_max:dx_max);
        for q1=ind_vecy
            %perform row shift
            dy=dy+1;
            curr_im2_temp=circshift(curr_im2,[q1,0]);
            curr_im2_temp_crop=curr_im2_temp(rect1(2):rect1(2)+rect1(4),:);
            
            dx=0;
            for q2=ind_vecx
                %perform column shift
                dx=dx+1;
                curr_im2_temp2=circshift(curr_im2_temp_crop,[0,q2]);
                
                %clip out correct region and normalize
                sub_im2_raw=curr_im2_temp2(:,rect1(1):rect1(1)+rect1(3));
                sub_im2_mean=sub_im2_raw-mean(sub_im2_raw(:));
                sub_im2=sub_im2_mean/sqrt(sum(sub_im2_mean(:).^2));
                
                %calculate correlation coefficient
                cc_mat1(dy,dx)=sum(sub_im1(:).*sub_im2(:));
            end
        end
        
        %find shift values
        [y0,x0]=find(cc_mat1==max(cc_mat1(:)),1,'first');
        
        %check if on edge
        if or(or(y0==1,y0==2*dx_max+1),or(x0==1,x0==2*dx_max+1))
            disp(['Warning: shift limit reached on frame ' num2str(i) '.'])
            bnd_warn=1;
        else
            bnd_warn=0;
        end
        
        %break loop the translation gets too high
        if dx_max>(2*dx_max0)
            bnd_warn=0;
            disp(['Warning: adaptive shift limit reached on frame ' num2str(i) '.'])
            disp(['Could not find a suitable shift value for this frame, all future'])
            disp(['frame shifts may be compromised.'])
        end
    end
    
    %record final values  
    X_shift(i)=ind_vecx(x0);
    Y_shift(i)=ind_vecy(y0);
end

set(handles.text_progress,'String',['Writing image files ...'])
drawnow
%write new image files
%for k=find(master_vec+slave_vec)
for k=find(types_vec~=0)'
    curr_file=dir_base(k).name;
    [~,name_temp,~]=fileparts(curr_file);
    basename=[name_temp '_aligned.tif'];
    
    for i=1:N
        %perform shift of current images
        shift_im=circshift(imread(curr_file,i),[Y_shift(i),X_shift(i)]);
        
        %save image
        image_save(shift_im,basename)
    end
end
set(handles.text_progress,'String',['Finished.'])
            

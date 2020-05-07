function varargout = stack_maker(varargin)
% STACK_MAKER MATLAB code for stack_maker.fig
%      STACK_MAKER, by itself, creates a new STACK_MAKER or raises the existing
%      singleton*.
%
%      H = STACK_MAKER returns the handle to a new STACK_MAKER or the handle to
%      the existing singleton*.
%
%      STACK_MAKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_MAKER.M with the given input arguments.
%
%      STACK_MAKER('Property','Value',...) creates a new STACK_MAKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_maker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_maker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_maker

% Last Modified by GUIDE v2.5 02-Sep-2015 02:39:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_maker_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_maker_OutputFcn, ...
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


% --- Executes just before stack_maker is made visible.
function stack_maker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_maker (see VARARGIN)

handles.output = hObject;

global folder_names
global file_names
global dirname

folder_names={};
file_names={};
dirname={};

%default motifs
set(handles.edit_submotif,'String','Pos*')
set(handles.edit_imagemotif,'String','img*.tif')

checkbox_subdirectories_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = stack_maker_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_directory.
function pushbutton_directory_Callback(hObject, eventdata, handles)
global folder_names
global file_names
global dirname

%get basal directory
dirname=uigetdir(pwd,'Select Root Image Directory');
if dirname==0
    dirname={};
    return
end
cd(dirname)

%reset text fields
set(handles.text_founddirs,'String',' ')
set(handles.text_foundfiles,'String',' ')
set(handles.text_progress_folder,'String',' ')
set(handles.text_progress_file,'String',' ')
set(handles.text_progress_type,'String',' ')

edit_submotif_Callback(hObject, eventdata, handles)
edit_imagemotif_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_convert.
function pushbutton_convert_Callback(hObject, eventdata, handles)
global folder_names
global file_names
global dirname

if isempty(dirname)
    warndlg('Select a root directory in which images files or sub-directories are located before converting.','Select Directory First')
    return
end

edit_submotif_Callback(hObject, eventdata, handles)
edit_imagemotif_Callback(hObject, eventdata, handles)

if isempty(folder_names)
    convert_images(hObject, eventdata, handles, file_names,folder_names,dirname)
else
    %handle each sub-directory
    for i=1:length(folder_names)
        
        %parse file names
        file_names_temp=dir(fullfile(dirname,fullfile(folder_names(i).name,get(handles.edit_imagemotif,'String'))));
        
        convert_images(hObject, eventdata, handles, file_names_temp,folder_names(i).name,dirname)
    end
end


function convert_images(hObject, eventdata, handles, file_names, folder_name,dirname)

%get the number and names of channels
name_struct=parse_input(get(handles.edit_imagemotif,'String'),'*');

%characters to remove from name
rm_text={' ',',','_','0','1','2','3','4','5','6','7','8','9',name_struct{:}};

for j=1:length(file_names)
    %get file name
    remain=file_names(j).name;
    
    %remove characters
    for i=1:length(rm_text)
        r_temp = strfind(remain, rm_text{i});
        
        %create filler
        filler=[];
        for k=1:length(rm_text{i})
            filler=[filler '#'];
        end
        
        %introduce filler
        for k=1:length(r_temp)
            remain(r_temp(k):r_temp(k)+length(rm_text{i})-1)=filler;
        end
    end
    %remove filler
    file_names(j).channel=[[name_struct{1:end-1}] '_' remain(remain~='#')];
end

%find number of channels
n_channels=length(unique({file_names(:).channel}));
dt=round(n_channels/length(file_names)*(file_names(end).datenum-file_names(1).datenum)*86400*10)/10;

if get(handles.checkbox_framecrop,'Value')
    start_fr=str2num(get(handles.edit_startframe,'String'));
    end_fr=str2num(get(handles.edit_endframe,'String'));
else
    start_fr=1;
    end_fr=length(file_names);
end

%check for RGB
rgb_check=1;

%for j=1:length(file_names)
for j=start_fr:end_fr
    %get channel type
    fname=file_names(j).channel;
    
    %construct output name
    if isempty(folder_name)
        %construct output name
        out_name=fullfile(dirname,[fname '_' num2str(dt) 's.tif']);
        %get input name
        in_name=fullfile(dirname,file_names(j).name);
    else
        %construct output name
        out_name=fullfile(dirname,[folder_name '_' fname '_' num2str(dt) 's.tif']);
        %get input name
        in_name=fullfile(dirname,fullfile(folder_name,file_names(j).name));
    end
    
    %write output image
    try
        %check to make sure images are not RGB
        if rgb_check
            info1=imfinfo(in_name);
            if ~strcmp(info1(1).ColorType,'grayscale')
            
                temp1=questdlg('These images seem to be non-grayscale (e.g. RGB).  Do you want to convert them to grayscale?', ...
                         'Convert Images', ...
                         'Yes', 'No','Yes');
                
                if strcmp(temp1,'Yes')
                    rgb_convert=1;
                else
                    rgb_convert=0;
                end    
            end
            rgb_check=0;
        end
        
        %load in image
        im_temp=imread(in_name);
        
        if rgb_convert
            im_temp2=zeros(size(im_temp,1),size(im_temp,2));
            for q1=1:size(im_temp,3)
                im_temp2=im_temp2+double(im_temp(:,:,q1)).^2;
            end
            im_temp2=uint16(sqrt(1/3*im_temp2));
        else
            im_temp2=im_temp;
        end
        
        if get(handles.checkbox_invert,'Value')
            image_save(imcomplement(im_temp2),out_name)
        else
            image_save(im_temp2,out_name)
        end
        
        if get(handles.checkbox_delete,'Value')
            delete(in_name)
        end
        
    catch er1
        info1=imfinfo(out_name);
        info2=imfinfo(in_name);
        if (info1(1).FileSize + info2(1).FileSize)>(2^32-1)
            disp('Oh jeez ... so ... MATLAB, stupidly, limits TIF file sizes to ~4GB, and you''ve now')
            disp('exceeded that limit. Time to use ImageJ / Fiji to turn these images into a stack.')
            disp('The ~4GB image stack that was just created is still present though.')
            break
        else
            disp(['Was unable to save image ' in_name ' to the stack.'])
            disp('File may be in use, or incorrectly formatted.')
            disp(' ')
            disp(er1.message)
            continue
        end
    end
    
    if isempty(folder_name)
        set(handles.text_progress_folder,'String',' ')
    else
        set(handles.text_progress_folder,'String',['Folder: ' folder_name])
    end
    set(handles.text_progress_file,'String',['File: ' num2str(j)]);
    set(handles.text_progress_type,'String',['Type: ' fname]);
    drawnow
end
fclose all;
set(handles.text_progress_folder,'String',['Finished, with mean frame time(s): ' num2str(dt)])


% --- Executes on button press in checkbox_subdirectories.
function checkbox_subdirectories_Callback(hObject, eventdata, handles)
global folder_names
if get(handles.checkbox_subdirectories,'Value')
    set(handles.text_subtext,'enable','on')
    set(handles.edit_submotif,'enable','on')
    set(handles.text_foundfiles,'String',' ')
    edit_submotif_Callback(hObject, eventdata, handles)
else
    set(handles.text_subtext,'enable','off')
    set(handles.edit_submotif,'enable','off')
    set(handles.text_founddirs,'String',' ')
    folder_names={};
end
edit_imagemotif_Callback(hObject, eventdata, handles)


function edit_submotif_Callback(hObject, eventdata, handles)
global dirname
global folder_names

if isempty(get(handles.edit_submotif,'String'))
    set(handles.edit_submotif,'String','*')
end

if and(~isempty(dirname),get(handles.checkbox_subdirectories,'Value'))
    %find files
    folder_names_temp=dir(fullfile(dirname,[get(handles.edit_submotif,'String')]));
    folder_names={};

    if ~isempty(folder_names_temp)
        q=0;
        for i=1:length(folder_names_temp)
            if and(folder_names_temp(i).isdir,~isempty(folder_names_temp(i).name))
                q=q+1;
                folder_names(q).name=folder_names_temp(i).name;
            end
        end
                
        set(handles.text_founddirs,'String',['found ' num2str(length(folder_names)) ' sub-directories'])
    else
        set(handles.text_founddirs,'String',['found 0 sub-directories'])
    end
end

% --- Executes during object creation, after setting all properties.
function edit_submotif_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_imagemotif_Callback(hObject, eventdata, handles)
global dirname
global file_names

if isempty(get(handles.edit_imagemotif,'String'))
    set(handles.edit_imagemotif,'String','*')
end

if ~isempty(dirname)
    %find files
    file_names_temp=dir(fullfile(dirname,[get(handles.edit_imagemotif,'String')]));
    file_names={};
    
    if ~isempty(file_names_temp)
        q=0;
        for i=1:length(file_names_temp)
            if and(~file_names_temp(i).isdir,~isempty(file_names_temp(i).name))
                q=q+1;
                file_names(q).name=file_names_temp(i).name;
                file_names(q).datenum=file_names_temp(i).datenum;
            end
        end
                
        if get(handles.checkbox_subdirectories,'value')
            set(handles.text_foundfiles,'String',['found ' num2str(length(file_names)) ' image files, in root directory'])
        else
            set(handles.text_foundfiles,'String',['found ' num2str(length(file_names)) ' image files'])
        end
    else
        if get(handles.checkbox_subdirectories,'value')
            set(handles.text_foundfiles,'String',['found 0 image files'])
        else
            set(handles.text_foundfiles,'String',['found 0 image files, in root directory'])
        end
    end
end

% --- Executes during object creation, after setting all properties.
function edit_imagemotif_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_delete.
function checkbox_delete_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_invert.
function checkbox_invert_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_framecrop.
function checkbox_framecrop_Callback(hObject, eventdata, handles)
if ~get(handles.checkbox_framecrop,'Value')
    set(handles.edit_startframe,'enable','off')
    set(handles.edit_endframe,'enable','off')
    set(handles.text16,'enable','off')
    set(handles.text17,'enable','off')
else
    set(handles.edit_startframe,'enable','on')
    set(handles.edit_endframe,'enable','on')
    set(handles.text16,'enable','on')
    set(handles.text17,'enable','on')
end


function edit_startframe_Callback(hObject, eventdata, handles)
global file_names
if str2num(get(handles.edit_startframe,'String'))<1
    set(handles.edit_startframe,'String',num2str(1))
end
if ~isempty(file_names)
    if str2num(get(handles.edit_startframe,'String'))>length(file_names)
        set(handles.edit_startframe,'String',num2str(length(file_names)-1))
        set(handles.edit_endframe,'String',num2str(length(file_names)))
    end
else
    set(handles.edit_startframe,'String','1')
end
edit_endframe_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_startframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_endframe_Callback(hObject, eventdata, handles)
global file_names

if str2num(get(handles.edit_endframe,'String'))>length(file_names)
    set(handles.edit_endframe,'String',num2str(length(file_names)))
end
if str2num(get(handles.edit_endframe,'String'))<str2num(get(handles.edit_startframe,'String'))
    set(handles.edit_endframe,'String',get(handles.edit_startframe,'String'))
end
% --- Executes during object creation, after setting all properties.
function edit_endframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%{
% --- Executes on button press in pushbutton_stacklength.
function pushbutton_stacklength_Callback(hObject, eventdata, handles)
global dirname

%get the basal tif stack name.
if get(handles.checkbox_subdirectories,'Value')
    base_motif=[get(handles.edit_submotif,'String') get(handles.edit_imagemotif,'String')];
else
    base_motif=[get(handles.edit_imagemotif,'String')];
end

%get file names
if isempty(dirname)
    dir0=dir(base_motif);
else
    dir0=dir(fullfile(dirname,base_motif));
end

%get frame numbers for re-save
start_fr=str2num(get(handles.edit_startframe,'String'));
end_fr=str2num(get(handles.edit_endframe,'String'));

%process files
q=0;
for n=1:length(dir0)
    if ~dir0(n).isdir
        q=q+1;
        
        set(handles.text_progress_folder,'String',['File: ' dir0(n).name]);
        drawnow
        
        %loop over images
        for i=start_fr:end_fr
            set(handles.text_progress_type,'String',['Frame: ' num2str(i)]);
            drawnow
            
            %get image
            temp1=imread(dir0(n).name,i);
            
            image_save(temp1,['stack_maker_temp.tif'])
        end
        
        %delete old file
        delete(dir0(n).name)
        
        %rename new file
        er1=0;
        f=0;
        while and(er1==0,f<11)
            try
                movefile(['stack_maker_temp.tif'],dir0(n).name)
                er1=1;
            catch
                f=f+1;
                pause(1*rand)
            end
        end
    end
end

set(handles.text_progress_folder,'String',['Finished cropping ' num2str(q) ' files.'])
set(handles.text_progress_file,'String',[' '])
set(handles.text_progress_type,'String',[' '])
fclose all;
%}

function varargout = contour_parameters(varargin)
% CONTOUR_PARAMETERS MATLAB code for contour_parameters.fig
%      CONTOUR_PARAMETERS, by itself, creates a new CONTOUR_PARAMETERS or raises the existing
%      singleton*.
%
%      H = CONTOUR_PARAMETERS returns the handle to a new CONTOUR_PARAMETERS or the handle to
%      the existing singleton*.
%
%      CONTOUR_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTOUR_PARAMETERS.M with the given input arguments.
%
%      CONTOUR_PARAMETERS('Property','Value',...) creates a new CONTOUR_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before contour_parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to contour_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help contour_parameters

% Last Modified by GUIDE v2.5 27-Jul-2015 19:24:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contour_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @contour_parameters_OutputFcn, ...
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


% --- Executes just before contour_parameters is made visible.
function contour_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contour_parameters (see VARARGIN)

% Choose default command line output for contour_parameters
handles.output = hObject;

global dl
global res
global offset

set(handles.edit_linewidth,'String',num2str(dl))
set(handles.edit_res,'String',num2str(res))
set(handles.edit_offset,'String',num2str(offset))

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contour_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = contour_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_linewidth_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_linewidth,'String'))<1
    set(handles.edit_linewidth,'String','1')
end
% --- Executes during object creation, after setting all properties.
function edit_linewidth_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_res_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_res,'String'))<1
    set(handles.edit_res,'String','1')
end
% --- Executes during object creation, after setting all properties.
function edit_res_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_offset_Callback(hObject, eventdata, handles)
temp1=str2num(get(handles.edit_offset,'String'));
if or(temp1<=0,temp1>1)
    set(handles.edit_offset,'String','0.5')
end
% --- Executes during object creation, after setting all properties.
function edit_offset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
global dl
global res
global offset

dl=str2num(get(handles.edit_linewidth,'String'));
res=str2num(get(handles.edit_res,'String'));
offset=str2num(get(handles.edit_offset,'String'));

close contour_parameters

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
close contour_parameters

function varargout = noise_filt_parameters(varargin)
% NOISE_FILT_PARAMETERS MATLAB code for noise_filt_parameters.fig
%      NOISE_FILT_PARAMETERS, by itself, creates a new NOISE_FILT_PARAMETERS or raises the existing
%      singleton*.
%
%      H = NOISE_FILT_PARAMETERS returns the handle to a new NOISE_FILT_PARAMETERS or the handle to
%      the existing singleton*.
%
%      NOISE_FILT_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NOISE_FILT_PARAMETERS.M with the given input arguments.
%
%      NOISE_FILT_PARAMETERS('Property','Value',...) creates a new NOISE_FILT_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before noise_filt_parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to noise_filt_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help noise_filt_parameters

% Last Modified by GUIDE v2.5 27-Jul-2015 19:21:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @noise_filt_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @noise_filt_parameters_OutputFcn, ...
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


% --- Executes just before noise_filt_parameters is made visible.
function noise_filt_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to noise_filt_parameters (see VARARGIN)

% Choose default command line output for noise_filt_parameters
handles.output = hObject;

global rel_sz
global rel_sigma

set(handles.edit_sz,'String',num2str(rel_sz))
set(handles.edit_sigma,'String',num2str(rel_sigma))

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes noise_filt_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = noise_filt_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_sz_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_sz,'String'))<3
    set(handles.edit_sz,'String','3')
end
% --- Executes during object creation, after setting all properties.
function edit_sz_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_sigma_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_sigma,'String'))<0.1
    set(handles.edit_sigma,'String','0.1')
end
% --- Executes during object creation, after setting all properties.
function edit_sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
global rel_sz
global rel_sigma

rel_sz=str2num(get(handles.edit_sz,'String'));
rel_sigma=str2num(get(handles.edit_sigma,'String'));

close noise_filt_parameters

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
close noise_filt_parameters

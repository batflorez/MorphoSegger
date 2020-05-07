function varargout = morphometrics(varargin)
%
% morphometrics v2.0
% Updated version by Andres Florez - 04/20/20
% **************************************************************************
%morphometrics v0.98
% (c) Tristan Ursell 2012,2013,2014.
%
% Please see the 'help' button for specific guidance.
%
%Acknowledgements:
% Carolina Tropini for help with the initial GUI.
% Timothy Lee for significant upgrades in the speed of code execution.
% Continued support from the Huang Group at Stanford
% (whatislife.stanford.edu)



% Last Modified by GUIDE v2.5 07-Oct-2016 19:42:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @morphometrics_OpeningFcn, ...
    'gui_OutputFcn',  @morphometrics_OutputFcn, ...
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

function morphometrics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to morphometrics (see VARARGIN)

fclose all;

%check for image processing toolbox
if ~license('test','image_toolbox')
    error('This version of MATLAB does not have the Image Processing Toolbox installed; Morphometrics cannot run without it. Sorry.')
end

%check matlab version
temp1=version;
ind1=strfind(temp1,'R');
if str2num(temp1(ind1+1:ind1+4))<2012
    error('Morphometrics requires Matlab release 2012b or later.')
end

%handle preferences (from file at first, then defaults if no file)
%Defaults based on wt E coli at 80 nm/pixel (4-2-2012)
global home_path %return path for common file locations
[home_path,~,~]=fileparts(mfilename('fullpath'));

%ask to add morphmetrics home folder to matlab path
try
    fid1=fopen(fullfile(home_path,'ask_add_to_path.txt'));
    temp1=textscan(fid1,'%u');
    fclose(fid1);
catch
    %turn off path add question dialog
    fid1= fopen(fullfile(home_path,'ask_add_to_path.txt'),'w+');
    fprintf(fid1,'%u',1);
    fclose(fid1);
    
    fid1=fopen(fullfile(home_path,'ask_add_to_path.txt'));
    temp1=textscan(fid1,'%u');
    fclose(fid1);
end

if temp1{1}
    ButtonName = questdlg('Add the Morphometrics home folder to the Matlab path?', ...
        'Add Morphometrics to Path', ...
        'Yes', 'No', 'Do not ask again', 'Yes');
    
    switch ButtonName,
        case 'Yes',
            menu_addpath_Callback(hObject, eventdata, handles)
        case 'No',
            warndlg('When everything is on fire, your teeth are gnashing, and the red velvet cupcake you just paid $9 for at some boutique bakery is in a nasty puddle -- you''ll know why.', 'Suit yourself.');
        case 'Do not ask again',
            %turn off path add question dialog
            fid= fopen(fullfile(home_path,'ask_add_to_path.txt'),'w+');
            fprintf(fid,'%u',0);
            fclose(fid);
    end
end

%version control check (Not necessary since the software is not going to be
%updated anymore - added by Andres Florez on 04/16/2020
% try
%     fid2=fopen(fullfile(home_path,'update_on_start.txt'));
%     temp2=textscan(fid2,'%u');
%     fclose(fid2);
%     if temp2{1}
%         menu_update_Callback(hObject, eventdata, handles)
%     end
% catch
%     menu_update_Callback(hObject, eventdata, handles)
% end

enable_buttons(handles)

%initialize and load all global variables and parameters
global stopV
stopV=0;  %set stop indicator

global testV
testV=0;  %set test indicator

global f_canny_th1 %lower canny threshold, triggers edge start
global f_canny_th2 %upper canny threshold, triggers edge end
global f_canny_sig %canny sigma factor, determines smoothness

global f_resize %image resizing factor
global f_back %background removal filter size
global f_histlim %histogram percentage removal for contrasting
global f_areamin %minimum object area
global f_areamax %maximum object area
global f_hmin %primary watershed threshold
global f_hmin_split %secondary pixel distance threshold
global f_gstd %initial Gaussian smoothing filter STD
global f_r_int %Gaussing force smoothing STD
global f_pert_same %percent overlap for cells to be related
global f_frame_diff %number of frames in back/front to look for related cells

global f_relsz %relative noise filter disk size (diameter)
global f_int_rej %false positive rejection parameter
%global f_seg_min %minimum relative variacne required for any segmentation
global f_seg_min %histogram tail weight used to determine if any segmentation happens
global f_contmin %minimum number of points a contour may have
global f_acute %parameter that controls removal of acute angles in contours
global f_pixelprox %number of pixels below which two cells are considered proximal
global f_curve_std %standard deviation (in pixels) of guassian convolution curvature filter

global f_rect %coordinates of the background ROI rectangle for adaptive thresholding

global f_Nit %number of contour fitting iterations
global f_kint %rate of contour fitting

global v_method %sets the method variable (gradient, laplacian, adaptive threshold)
global v_imtype %sets type of image being used (phase dark, interior fluor, peripheral fluor)

%it seems that force smoothing doesn't do anything useful, thus turning off
%this parameter
set(handles.edit_force_smooth,'Visible','off')
set(handles.text_force_smooth,'Visible','off')

try
    load(fullfile(home_path,'Morphometrics_prefs.mat'),'f_*','v_*')
    
    set(handles.edit_cannylow,'String',num2str(f_canny_th1))
    set(handles.edit_cannyhigh,'String',num2str(f_canny_th2))
    set(handles.edit_cannysigma,'String',num2str(f_canny_sig))
    
    set(handles.edit_resize,'String',num2str(f_resize))
    set(handles.edit_back,'String',num2str(f_back))
    set(handles.edit_histlim,'String',num2str(f_histlim))
    set(handles.edit_min_area,'String',num2str(f_areamin))
    set(handles.edit_max_area,'String',num2str(f_areamax))
    set(handles.edit_smooth,'String',num2str(f_gstd))
    set(handles.edit_force_smooth,'String',num2str(f_r_int))
    set(handles.edit_overlap,'String',num2str(f_pert_same))
    set(handles.edit_frame,'String',num2str(f_frame_diff))
    
    switch v_method
        case 1
            checkbox_seg_grad_Callback(hObject, eventdata, handles)
        case 2
            checkbox_seg_lap_Callback(hObject, eventdata, handles)
        case 3
            checkbox_seg_thres_Callback(hObject, eventdata, handles)
        case 4
            checkbox_seg_canny_Callback(hObject, eventdata, handles)
    end
    
    switch v_imtype
        case 1
            set(handles.checkbox_phase,'Value',1);
            set(handles.checkbox_fluor_int,'Value',0);
            set(handles.checkbox_fluor_periphery,'Value',0);
        case 2
            set(handles.checkbox_phase,'Value',0);
            set(handles.checkbox_fluor_int,'Value',1);
            set(handles.checkbox_fluor_periphery,'Value',0);
        case 3
            set(handles.checkbox_phase,'Value',0);
            set(handles.checkbox_fluor_int,'Value',0);
            set(handles.checkbox_fluor_periphery,'Value',1);
    end
    
    set(handles.checkbox_indpt,'Value',v_indpt)
    set(handles.checkbox_prox,'Value',v_prox)
    set(handles.checkbox_save,'Value',v_save)
    set(handles.checkbox_falsepos,'Value',v_falsepos)
    set(handles.checkbox_advanced,'Value',v_advanced)
    set(handles.checkbox_seed,'Value',v_seed)
    set(handles.checkbox_manualtrack,'Value',v_manualtrack)
    set(handles.checkbox_mip_persist,'Value',v_persist)
    set(handles.checkbox_simplethres,'Value',v_simplethres)
    set(handles.checkbox_colorcode,'Value',v_colorcode)
    set(handles.checkbox_excludeedge,'Value',v_exclude)
    
    set(handles.checkbox_keeptemp,'Value',0)
    set(handles.checkbox_savetest,'Value',0)
    
    set(handles.checkbox_meshmotif,'Value',v_meshmotif)
    set(handles.checkbox_mm_mesh,'Value',v_mm_mesh)
    set(handles.checkbox_mt_mesh,'Value',v_mt_mesh)
    
    set(handles.edit_primary,'String',num2str(f_hmin))
    set(handles.edit_secondary,'String',num2str(f_hmin_split))
catch
    set(handles.text_process,'String','The Morphometrics preference file seems to be missing or corrupted. A new one was created with default values.')
    
    %***************
    %DEFAULT VALUES*
    f_canny_th1=0.4; %lower canny threshold, triggers edge start
    f_canny_th2=0.8; %upper canny threshold, triggers edge end
    f_canny_sig=1; %canny sigma factor, determines smoothness
    
    f_resize=1; %image resizing factor
    f_back=0; %background removal filter size
    f_histlim=0.001; %histogram percentage removal for constrasting
    f_areamin=200; %minimum object area
    f_areamax=20000; %maximum object area
    f_hmin=0.1; %primary threshold (for gradient mode)
    f_hmin_split=5; %secondary splitting based on pixel distance transform
    f_gstd=0; %Gaussian image smoothing STD
    f_r_int=1; %Guassian force smoothing STD
    f_pert_same=0.55; %percent overlap for cells to be related
    f_frame_diff=4; %number of frames in back/front to look for related cells
    
    f_relsz=3; %relative noise filter disk size (diameter)
    f_int_rej=1; %false positive rejection parameter
    f_seg_min=0.02; %histogram tail weight used to determine if any segmentation happens
    f_contmin=10; %minimum number of points a contour may have
    f_acute=1; %parameter that controls removal of acute angles in contours
    f_pixelprox=5; %number of pixels below which two cells are considered proximal
    f_curve_std=2; %standard deviation (in pixels) of guassian convolution curvature filter
    
    f_rect=0; %coordinates of the background ROI rectangle for adaptive thresholding
    
    v_method = 1; %set segmentation method (1=gradient, 2=laplacian, 3=adaptive threshold)
    v_imtype = 1; %set imgage type method (1=phase dark, 2=intertio fluor, 3= peripheral fluor)
    
    v_indpt = 0; %sets independent images analysis checkbox value
    v_prox = 0; %sets cell proximity checkbox value
    v_save = 1; %sets save checkbox value
    v_falsepos = 1; %sets reject false positive checkbox
    v_advanced = 0; %choose between advanced and basic mode
    v_seed = 0; %sets whether contours are seeded by previous frame
    v_manualtrack = 0; %sets whether frame to frame tracking is manually curated
    v_persist = 0; %sets whether to use persistent background ROI during adaptive thresholding
    v_simplethres = 0; %sets whether to use a simple intensity threshold
    v_colorcode = 0; %sets whether the output segmentation image is colored coded by rejection type
    v_exclude = 1; %set whether objects on the image edge are excluded from segmentation
    v_nocontours=1; %set whether objects will have contours (1) or just segmented pixels (0)
    
    v_meshmotif = 0; %sets mesh motif checkbox value
    v_mm_mesh = 1; %sets meshing to 'morphometrics' mesh
    v_mt_mesh = 0; %sets meshing to 'microbe_tracker' mesh
    
    set(handles.checkbox_keeptemp,'Value',0);
    set(handles.checkbox_savetest,'Value',0);
    
    %Don't change these values unless you know what you are doing
    f_Nit=300; %number of fitting iterations
    f_kint=0.25; %rate of contour fitting
    %***************
    
    set(handles.edit_cannylow,'String',num2str(f_canny_th1))
    set(handles.edit_cannyhigh,'String',num2str(f_canny_th2))
    set(handles.edit_cannysigma,'String',num2str(f_canny_sig))
    
    set(handles.edit_resize,'String',num2str(f_resize))
    set(handles.edit_back,'String',num2str(f_back))
    set(handles.edit_histlim,'String',num2str(f_histlim))
    set(handles.edit_min_area,'String',num2str(f_areamin))
    set(handles.edit_max_area,'String',num2str(f_areamax))
    set(handles.edit_primary,'String',num2str(f_hmin))
    set(handles.edit_secondary,'String',num2str(f_hmin_split))
    set(handles.edit_smooth,'String',num2str(f_gstd))
    set(handles.edit_force_smooth,'String',num2str(f_r_int))
    set(handles.edit_overlap,'String',num2str(f_pert_same))
    set(handles.edit_frame,'String',num2str(f_frame_diff))
    
    switch v_method
        case 1
            checkbox_seg_grad_Callback(hObject, eventdata, handles)
        case 2
            checkbox_seg_lap_Callback(hObject, eventdata, handles)
        case 3
            checkbox_seg_thres_Callback(hObject, eventdata, handles)
    end
    
    switch v_imtype
        case 1
            set(handles.checkbox_phase,'Value',1);
            set(handles.checkbox_fluor_int,'Value',0);
            set(handles.checkbox_fluor_periphery,'Value',0);
        case 2
            set(handles.checkbox_phase,'Value',0);
            set(handles.checkbox_fluor_int,'Value',1);
            set(handles.checkbox_fluor_periphery,'Value',0);
        case 3
            set(handles.checkbox_phase,'Value',0);
            set(handles.checkbox_fluor_int,'Value',0);
            set(handles.checkbox_fluor_periphery,'Value',1);
    end
    
    set(handles.checkbox_indpt,'Value',v_indpt);
    set(handles.checkbox_prox,'Value',v_prox);
    set(handles.checkbox_save,'Value',v_save);
    set(handles.checkbox_falsepos,'Value',v_falsepos);
    set(handles.checkbox_keeptemp,'Value',0);
    set(handles.checkbox_advanced,'Value',v_advanced);
    set(handles.checkbox_seed,'Value',v_seed);
    set(handles.checkbox_manualtrack,'Value',v_manualtrack);
    set(handles.checkbox_mip_persist,'Value',v_persist);
    set(handles.checkbox_simplethres,'Value',v_simplethres);
    set(handles.checkbox_colorcode,'Value',v_colorcode);
    set(handles.checkbox_excludeedge,'Value',v_exclude);
    set(handles.checkbox_nocontours,'Value',v_nocontours);

   
    set(handles.checkbox_meshmotif,'Value',v_meshmotif);
    set(handles.checkbox_mm_mesh,'Value',v_mm_mesh);
    set(handles.checkbox_mt_mesh,'Value',v_mt_mesh);
    
    %create default prefs file
    save(fullfile(home_path,'Morphometrics_prefs.mat'),...
        'f_canny_th1',...
        'f_canny_th2',...
        'f_canny_sig',...
        'f_resize',...
        'f_back',...
        'f_histlim',...
        'f_areamin',...
        'f_areamax',...
        'f_hmin',...
        'f_hmin_split',...
        'f_gstd',...
        'f_r_int',...
        'f_pert_same',...
        'f_frame_diff',...
        'f_Nit',...
        'f_kint',...
        'f_relsz',...
        'f_int_rej',...
        'f_seg_min',...
        'f_contmin',...
        'f_acute',...
        'f_pixelprox',...
        'f_curve_std',...
        'f_rect',...
        'v_method',...
        'v_imtype',...
        'v_indpt',...
        'v_prox',...
        'v_save',...
        'v_falsepos',...
        'v_advanced',...
        'v_seed',...
        'v_manualtrack',...
        'v_persist',...
        'v_simplethres',...
        'v_colorcode',...
        'v_exclude',...
        'v_nocontours',...
        'v_meshmotif',...
        'v_mm_mesh',...
        'v_mt_mesh');
    set(handles.text_process,'String','Parameter values saved to: Morphometrics_prefs.mat');
end

% Choose default command line output for morphometrics
handles.output = hObject;
set(handles.text_static,'String','*.tif')  %<--- change this to change the input file default name
set(handles.text_motif,'String','*motif*.tif')  %<--- change this to change the default motif name

global path1
global fnames
global mesh_motif_text
path1=0;
fnames=[];
mesh_motif_text=['*motif*CONTOURS.mat'];

f_rect=0; %coordinates of the background ROI

imagesc(imread('morphometrics_v2_f.jpg'),'Parent',handles.axes1);  %<--- change this if you want a different background picture
%set(handles.axes1,'Visible','on')
set(handles.axes1,'XTick',[])
set(handles.axes1,'YTick',[])
set(handles.axes1,'Box','on')

checkbox_indpt_Callback(hObject, eventdata, handles)
checkbox_prox_Callback(hObject, eventdata, handles)

checkbox_advanced_Callback(hObject, eventdata, handles)
guidata(hObject,handles)


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
tic
guidata(hObject,handles);
global path1
global fnames
global testV
global f_rect

global X_roi
global Y_roi

global f_resize

global v_method
global v_imtype

global stopV
stopV=0;

if path1==0
    warndlg('First load an input file, then try processing data.', 'Load Data First');
    return
end

%disable appropriate buttons
disable_buttons(handles)

%runs for multiple files
for i=1:length(fnames)
    clear frame cell %added by KC on 11/17/2015
    
    %number of stack images
    fname=fullfile(path1,fnames(i).name);
    Im_info=imfinfo(fname);
    set(handles.text_bitdepth,'String',[num2str(Im_info(1).BitDepth) ' bit'])
    set(handles.text_posbox,'String',[num2str(Im_info(1).Height) ' H  x  ' num2str(Im_info(1).Width) ' W'])
    
    Nim=length(Im_info);
    
    %setup roi mask
    mask1=ones(Im_info(1).Height,Im_info(1).Width);
    if and(i,exist('X_roi','var'))
        if length(X_roi)>3
            mask1=roipoly(mask1,X_roi,Y_roi);
        end
    end
    
    %get file name parts
    [pathstr, name, ~] = fileparts(fname);
    basename=fullfile(pathstr,name);
    
    %get rid of temporary files if they exist
    dir1=dir([basename '_G*.tif']);
    if ~isempty(dir1)
        delete([basename '_Gparent.tif']);
        delete([basename,'_Gsegt.tif']);
    end
    
    %update status text
    if length(fnames)>1
        set(handles.text_static,'String',[strtrim(fnames(i).name) ' (' num2str(i) ' of ' num2str(length(fnames)) ')']);
    else
        set(handles.text_static,'String',[strtrim(fnames(i).name)]);
    end
    
    set(handles.text_process,'String',['Starting segmentation of: ' strtrim(fnames(i).name)])
    title('Starting segmentation ...')
    drawnow
    
    %compare input file size with OS limits
    info0=imfinfo(fname);
    
    cond1=and(info0(1).BitDepth==8,info0(1).FileSize>2e9);
    cond2=and(info0(1).BitDepth==16,info0(1).FileSize>4e9);

    if and(or(cond1,cond2),~testV)
        error_q= questdlg([['The input image is too large to be handled by your operating system. '] ...
            ['We will try to fix this by limiting the maximum number of objects per frame to 255. ']...
            ['If there is more than 255 objects per frame, you must reduce the input file size either by cropping or frame number reduction.']]...
            , 'Input Data Too Large', ...
                         'Continue with 255 limit', 'Stop and reduce input file', 'Continue with 255 limit');
        
        if strcmp(error_q,'Continue with 255 limit')
            img_error=1;
        else
            set(handles.text_process,'String',['Please reduce input file size by cropping or frame reduction using (e.g.) ImageJ / FIJI.'])
            return
        end
    else
        img_error=0;
    end
    
    %loop over individual images
    for j=1:Nim
        set(handles.text_extrema,'String',[num2str(Im_info(j).MinSampleValue) ' Min,  ' num2str(Im_info(j).MaxSampleValue) ' Max'])
        
        %check for stop
        global stopV
        if stopV
            stopV=0;
            
            %delete temporary files
            delete([basename '_Gparent.tif']);
            delete([basename,'_Gsegt.tif']);
            
            set(handles.text_process,'String',['Segmentation stopped.'])
            title(['Segmentation stopped.'])
            return
        end
        
        %if this is only a test run
        if testV
            %get test frame number
            t_frame=round(str2num(get(handles.edit_test_frame,'String')));
                      
            Im0=imread(fullfile(path1,fnames(i).name),t_frame);
            
            %save temp image if desired (can then be individually processed)
            if get(handles.checkbox_savetest,'Value')
                %check for old temp files
                temp_name=[basename '_morphometrics_postprocess_presegment.tif'];
                dir1=dir(temp_name);
                
                if ~isempty(dir1)
                    delete(temp_name)
                end
                
                %write test image as a temp file
                image_save(Im0,temp_name)
            end
        else
            Im0=imread(fullfile(path1,fnames(i).name),j);
        end
        
        %perform saturation detection
        if Im_info(1).BitDepth==8
            sat_test=or(Im0==0,Im0==(2^8-1));
        else
            sat_test=or(Im0==0,Im0==(2^16-1));
        end
        
        sat_props=regionprops(sat_test,'Area');
        
        if any([[sat_props.Area]>=9])
            set(handles.text_saturation,'String','Saturation detected!')
        end
        
        %***********************************************
        %  Pre-processing for adaptive threshold
        %***********************************************
        if all([j==1,v_method==3,~get(handles.checkbox_simplethres,'Value')])
            %does not work with flourescence periphery
            if v_imtype==3
                warndlg('Adaptive Threshold segmentation is incompatible with Fluorescence (peripheral) images.', ' ');
                return
            end
            
            if or(f_rect==0,~get(handles.checkbox_mip_persist,'Value'))
                %create maximum intensity projection (mip)
                mip1=zeros(size(Im0));
                for k=1:Nim
                    temp1=imread(fullfile(path1,fnames(i).name),k);
                    mip1=max(cat(3,mip1,temp1),[],3);
                    title(['Calculating maximum intensity projection ... ' num2str(k) ' / ' num2str(Nim)])
                    drawnow
                end
                
                %display mip
                f1b=figure('units','normalized','outerposition',[0 0 1 1]);
                imagesc(mip1)
                axis equal tight
                title('Select a boxed background region from maximum intensity projection of the stack.')
                colormap(gray)
                f_rect=round(getrect);
                f_rect(f_rect<1)=1;
                
                while (f_rect(3)*f_rect(4))<25
                    title('Select a boxed background region from maximum intensity projection of the stack, it must be larger than 25 pixels.')
                    f_rect=round(getrect);
                    f_rect(f_rect<1)=1;
                end
                
                close(f1b)
            end
        elseif or(v_method~=3,and(v_method==3,get(handles.checkbox_simplethres,'Value')))
            f_rect=0;
        end
        %*******************************
        %*******************************
        
        %*******************************
        %     SEGMENTATION
        %*******************************
        %segment the images
        [segdata,Imag_out,Iseg_out,Im1,Icolor]=simply_segment(Im0,mask1,j,f_rect,handles);
        
        %save output segmentation data to main data structure
        if ~isempty(segdata)
            if segdata.num_objs==0
                frame(j).num_objs=0;
                frame(j).object=[];
            else
                frame(j)=segdata;
            end
        elseif and(j==1,get(handles.checkbox_seed,'value'))
            frame(j).num_objs=frame(1).num_objs;
        else
            frame(j).num_objs=0;
            frame(j).object=[];
        end
               
        %save output to stack for later processing
        try
            if img_error==0
                image_save(Iseg_out,[basename '_Gparent.tif'],20)
                image_save(Imag_out,[basename '_Gsegt.tif'],20)
            else
                image_save(uint8(Iseg_out),[basename '_Gparent.tif'],20)
                image_save(uint8(255*mat2gray(Imag_out)),[basename '_Gsegt.tif'],20)
            end
        catch err1
            disp(err1.message)
            info1=imfinfo([basename '_Gparent.tif']);
            if and(img_error==0,info1(1).FileSize >4e9)
                disp(['Unfortunately, your operating system does not allow temporary files to be written that are greater'])
                disp(['than 4GB.  We suggest you note the current frame being processed, and crop your data set to'])
                disp(['that number of frames minus 1 and then re-run Morphometrics on that reduced data file.'])
                disp([' '])
                disp(['Sorry, this is not Morphometerics'' fault.'])
                disp([' '])
                disp(['However, we will try to fix this ... see the pop-up.'])
                disp([' '])
                img_error=1;
                
                if img_error
                    quest1 = questdlg('Is the number of objects in each frame less than 255?', ...
                        'Reduce Temporary Data', ...
                        'Yes (we can continue)', 'No (end of the line)', 'Yes (we can continue)');
                    
                    if strcmp(quest1,'Yes (we can continue)')
                        disp(['Great, we will now reformat the temporary data ... this may take a moment.'])
                        %resave as 8 bit instead of 16 bit
                        for h=1:j-1
                            temp_err_1=uint8(imread([basename '_Gparent.tif'],h));
                            temp_err_2=uint8(255*mat2gray(imread([basename '_Gsegt.tif'],h)));
                            
                            image_save(temp_err_1,[basename '_Gparent_err.tif'],20)
                            image_save(temp_err_2,[basename '_Gsegt_err.tif'],20)
                            
                            disp(['Converting ' num2str(h) ' of ' num2str(j) '.'])
                        end
                        
                        %delete old temporary files
                        delete([basename '_Gparent.tif']);
                        delete([basename '_Gsegt.tif']);
                        
                        %rename temporary temporary files
                        movefile([basename '_Gparent_err.tif'],[basename '_Gparent.tif'])
                        movefile([basename '_Gsegt_err.tif'],[basename '_Gsegt.tif'])
                    else
                        error(['Morphometrics has stopped -- please reduce the number of input frames.'])
                    end
                end
            elseif and(img_error==1,info1.FileSize >4e9)
                error(['Remember that temporary data file error from earlier? Even the reduced temporary data is too large.  Sorry, Morphometrics will stop now -- your OS is to blame.'])
            end
        end
        
        imagesc(Icolor)
        hold on
        for k=1:frame(j).num_objs
            if isfield(frame(j).object(k),'Xcent')
                text(frame(j).object(k).Xcent,frame(j).object(k).Ycent,num2str(k),'color',[0 0 1])
            end
        end
        axis equal tight
        colormap(gray)
        set(handles.axes1,'Box','on')
        hold off
        if all([j>1,testV==0,get(handles.checkbox_seed,'Value')])
            title(['Frame ' num2str(j) ' of ' num2str(Nim) ', calculating image gradients ...'])
            set(handles.text_process,'String',['Frame ' num2str(j) ' of ' num2str(Nim) ', calculating image gradients ...'])
        elseif testV
            title(['Frame ' num2str(t_frame) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects (finished preprocessing)'])
        else
            title(['Frame ' num2str(j) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects'])
            set(handles.text_process,'String',['Frame ' num2str(j) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects'])
        end
        drawnow
        
        clear Imag_out
        clear Iseg_out
        clear segdata
        %-------------------------------
        %-------------------------------
        if testV
            break
        end
    end
    
    %check to make sure something was found
    badset=0;
    if sum([frame.num_objs])<1
        set(handles.text_process,'String','No objects were found in this stack using current parameter values. Consider adjusting the parameter values.')
        testV=0;
        enable_buttons(handles)
        badset=1;
        
        %get rid of temporary files if they exist
        dir1=dir([basename '_G*.tif']);
        if ~isempty(dir1)
            delete([basename '_Gparent.tif']);
            delete([basename '_Gsegt.tif']);
        end
        
        break
    end
    
    if ~testV
        title(['Frame ' num2str(j) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects (finished preprocessing)'])
    end
    
    %*******************************
    %     LINEAGE
    %*******************************
    if all([~testV,Nim>1,~get(handles.checkbox_indpt,'Value'),~get(handles.checkbox_seed,'Value')])
        if ~get(handles.checkbox_manualtrack,'Value')
            %find unique cells in stack
            try
                global f_pert_same
                global f_frame_diff
                cell=family([basename '_Gparent.tif'],f_pert_same,0,f_frame_diff,handles);
            catch er1
                disp('An error was encountered while trying to create the cell lineage.')
                disp('Try adjusting parameter values in the ''Parameters'' pane, or try')
                disp('selecting ''Images are independent''.')
                disp(' ')
                disp('Analysis will continue using independent cell analysis on this stack.')
                disp(' ')
                disp(er1)
                
                set(handles.text_process,'String','Lineage error: Analysis will continue using independent cell analysis on this stack.')
                clear cell
            end
        else
            %MANUAL OBJECT TRACKING
            [cell,~]=manual_tracker(basename,fname,Nim,handles);
        end
    end
    
    %update frame data structure to contain cell lineage data
    if exist('cell','var')
        for j=1:length(cell)
            for k=1:length(cell(j).frames)
                %get frame number
                fr_num=cell(j).frames(k);
                
                %get bw_label number
                bw_num=cell(j).bw_label(k);
                
                frame(fr_num).object(bw_num).bw_label=bw_num;
                frame(fr_num).object(bw_num).cellID=j;
            end
        end
    end
    
    %update frame data structure to contain *seeded* cell lineage
    if get(handles.checkbox_seed,'Value')
        for j=1:Nim
            frame(j).num_objs=frame(1).num_objs;
            for k=1:frame(1).num_objs
                frame(j).object(k).bw_label=k;
                frame(j).object(k).cellID=k;
            end
        end
    end
    
    %give cells unique names when frames are independent
    if any([testV,Nim==1,get(handles.checkbox_indpt,'Value')])
        q=0;
        for j=1:Nim
            for k=1:frame(j).num_objs
                q=q+1;
                
                frame(j).object(k).bw_label=k;
                frame(j).object(k).cellID=q;
            end
            
            if testV
                break
            end
        end
    end
    
    %determine maximum cell number
    if ~testV
        Ncell=1;
        for j=1:Nim
            for k=1:frame(j).num_objs
                if isfield(frame(j).object(k),'cellID')
                    if frame(j).object(k).cellID>Ncell
                        Ncell=frame(j).object(k).cellID;
                    end
                end
            end
        end
    end
    
    %-------------------------------
    %-------------------------------
    
    %***********************************************
    % PROXIMITY TESTING & CONTOUR FITTING
    %***********************************************
    set(handles.text_process,'String',['Fitting contours ...']);
    drawnow
    
    if ~get(handles.checkbox_nocontours,'Value')
        for j=1:Nim
            if and(frame(j).num_objs~=0,isfield(frame(j).object,'cellID'))
                temp_frame=frame(j);
                if and(j>1,get(handles.checkbox_seed,'Value'))
                    temp_frame_m1=frame(j-1);
                    [frameout]=fit_contours(basename,temp_frame,j,Nim,handles,temp_frame_m1);
                else
                    [frameout]=fit_contours(basename,temp_frame,j,Nim,handles);
                end
                
                %check for stop
                global stopV
                if stopV
                    stopV=0;
                    
                    %delete temporary files
                    delete([basename '_Gparent.tif'])
                    delete([basename '_Gsegt.tif'])
                    enable_buttons(handles)
                    return
                end
                
                frame(j)=frameout;
                if testV
                    title(['Frame ' num2str(t_frame) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects (finished)'])
                else
                    title(['Frame ' num2str(j) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects (finished)'])
                end
            end
            
            if testV
                break
            end
        end
    end
    
    %output for saving test image
    if and(testV,get(handles.checkbox_savetest,'Value'))
        try
            time1=clock;
            testname=[basename '_morphometrics_test_' date '_' num2str(time1(4)) '_' num2str(time1(5)) '_' num2str(round(time1(6)))];
            imwrite(frame2im(getframe),[testname '.tif'],'compression','none')
            print(gcf,'-depsc2',[testname '.eps'])
        catch
            warndlg('Please move the Morphometrics window into your primary monitor and test again.','Window Off Screen')
            testV=0;
            enable_buttons(handles)
            return
        end
    end
    %-------------------------------
    %-------------------------------
  
    %get rid of temporary files if they exist
    dir1=dir([basename '_G*.tif']);
    if and(~get(handles.checkbox_keeptemp,'Value'),~isempty(dir1))
        delete([basename '_Gparent.tif']);
        delete([basename '_Gsegt.tif']);
    end
    
    if testV
        break
    end    
    %*******************************
    %     OUTPUT FILES
    %*******************************
    %assemble full data structure and write output data
    if and(~testV,get(handles.checkbox_save,'Value'))
        set(handles.text_process,'String','Saving data ...')
        guidata(hObject,handles)
        drawnow
        
        global f_canny_th1
        global f_canny_th2
        global f_canny_sig
        
        global f_resize
        global f_back
        global f_histlim
        global f_areamin
        global f_areamax
        global f_hmin
        global f_hmin_split
        global f_gstd
        global f_pert_same
        global f_frame_diff
        global f_r_int
        
        global f_relsz
        global f_int_rej
        global f_seg_min
        global f_contmin
        global f_acute
        global f_pixelprox
        global f_curve_std
        
        global f_rect
        
        global f_Nit
        global f_kint
        
        global v_method
        global v_imtype
        
        v_indpt=get(handles.checkbox_indpt,'Value');
        v_prox=get(handles.checkbox_prox,'Value');
        v_save=get(handles.checkbox_save,'Value');
        v_falsepos=get(handles.checkbox_falsepos,'Value');
        v_advanced=get(handles.checkbox_advanced,'Value');
        v_seed=get(handles.checkbox_seed,'Value');
        v_manualtrack=get(handles.checkbox_manualtrack,'Value');
        v_persist=get(handles.checkbox_mip_persist,'Value');
        v_simplethres=get(handles.checkbox_simplethres,'Value');
        v_colorcode=get(handles.checkbox_colorcode,'Value');
        v_exclude=get(handles.checkbox_excludeedge,'Value');
        v_nocontours=get(handles.checkbox_nocontours,'Value');
        
        v_meshmotif=get(handles.checkbox_meshmotif,'Value');
        v_mm_mesh=get(handles.checkbox_mm_mesh,'Value');
        v_mt_mesh=get(handles.checkbox_mt_mesh,'Value');
        
        outname=[basename '_' date '_CONTOURS.mat'];
        %clearvars -except frame cell outname v_* f_*
        %save(outname,'frame','cell','outname','v_*','f_*')
        save(outname,'frame','outname','Ncell','v_*','f_*')
        set(handles.text_process,'String','Finished.')
        set(handles.text_process,'String',['All relevant data has been saved to the file:' outname])
    else
        set(handles.text_process,'String',['Finished: ' fnames(i).name])
    end
end
%-------------------------------
%-------------------------------

if ~badset
    t1=toc;
    testV=0;
    set(handles.text_process,'String',['Finished in ' num2str(round(10*t1/60)/10) ' minutes.'])
    
    %enable appropriate buttons
    enable_buttons(handles)
end

guidata(hObject,handles)


% --- Executes on button press in open_file.
function open_file_Callback(hObject, eventdata, handles)
global home_path
global path1
global fnames

cd(home_path)

if path1==0
    [file0,path0]=uigetfile({'*.tif;*.tiff','TIF Image Files'},'Select Image Stack',home_path);
else
    [file0,path0]=uigetfile({'*.tif;*.tiff','TIF Image Files'},'Select Image Stack',path1);
end

if file0==0
    return
else
    path1=path0;
    fnames=[];
    fnames(1).name=file0;
end

%erase temporary files
delete(fullfile(home_path,'*_Gparent*'));
delete(fullfile(home_path,'*_Gsegt*'));

%get stack size
Im_info=imfinfo(fullfile(path1,fnames(1).name));
set(handles.text_bitdepth,'String',[num2str(Im_info(1).BitDepth) ' bit'])
set(handles.text_posbox,'String',[num2str(Im_info(1).Height) ' H  x  ' num2str(Im_info(1).Width) ' W'])
try
    set(handles.text_extrema,'String',[num2str(Im_info(1).MinSampleValue) ' Min,  ' num2str(Im_info(1).MaxSampleValue) ' Max'])
catch
    set(handles.text_extrema,'String',['intensity data not available'])
end
set(handles.text_saturation,'String',' ')

set(handles.text_process,'String','Ready')
set(handles.text_of,'String',['of ' num2str(length(Im_info))]);
set(handles.edit_test_frame,'String',num2str(1));
set(handles.text_static,'String',strtrim(fnames(1).name));

hold off
imagesc(imread(fullfile(path1,fnames(1).name),1))
colormap(gray)
axis equal tight

if ~strcmp(Im_info(1).ColorType,'grayscale')
    warndlg('This seems to be a non-grayscale image, please convert the stack to an 8 or 16 bit unsigned grayscale.', 'Non-Grayscale Image Detected');
end

%delete any ROI objects
pushbutton_deleteroi_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


% --- Executes on button press in open_mult.
function open_mult_Callback(hObject, eventdata, handles)
global home_path
global path1
global fnames

path1=uigetdir(home_path,'Select Folder');
if path1==0
    return
end

%get new file names
motif = get(handles.text_motif,'String');
fnames=dir(fullfile(path1,motif));

if isempty(fnames)
    warndlg('No files in this directory match the search string.','No Files Found')
    path1=0;
    fnames=[];
    return
end

%erase temporary files
delete(fullfile(home_path,'*Gparent*'));
delete(fullfile(home_path,'*Gsegt*'));

%turn off persisent background ROI
if get(handles.checkbox_seg_thres,'Value')
    set(handles.checkbox_mip_persist,'Value',0)
end

%get stack size
Im_info=imfinfo(fullfile(path1,fnames(1).name));
set(handles.text_bitdepth,'String',[num2str(Im_info(1).BitDepth) ' bit'])
set(handles.text_posbox,'String',[num2str(Im_info(1).Height) ' H  x  ' num2str(Im_info(1).Width) ' W'])
set(handles.text_extrema,'String',[num2str(Im_info(1).MinSampleValue) ' Min,  ' num2str(Im_info(1).MaxSampleValue) ' Max'])

set(handles.text_process,'String','Ready')
set(handles.text_of,'String',['of ' num2str(length(Im_info))]);
set(handles.edit_test_frame,'String',num2str(1));
set(handles.text_static,'String',[fnames(1).name ' (1 of ' num2str(length(fnames)) ')']);

hold off
imagesc(imread(fullfile(path1,fnames(1).name),1))
colormap(gray)
axis equal tight

%delete any ROI objects
pushbutton_deleteroi_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


%function disable buttons during running or testing
function disable_buttons(handles)
%disable appropriate buttons
set(handles.open_file,'enable','off')
set(handles.open_mult,'enable','off')

set(handles.run,'enable','off')
set(handles.pushbutton_prefs,'enable','off')
set(handles.pushbutton_load_prefs,'enable','off')
set(handles.pushbutton_test,'enable','off')
set(handles.pushbutton_roi,'enable','off')
set(handles.pushbutton_deleteroi,'enable','off')

set(handles.pushbutton_mesh,'enable','off')


%function disable buttons during running or testing
function enable_buttons(handles)
%disable appropriate buttons
set(handles.open_file,'enable','on')
set(handles.open_mult,'enable','on')

set(handles.run,'enable','on')
set(handles.pushbutton_prefs,'enable','on')
set(handles.pushbutton_load_prefs,'enable','on')
set(handles.pushbutton_test,'enable','on')
set(handles.pushbutton_roi,'enable','on')
set(handles.pushbutton_deleteroi,'enable','on')

set(handles.pushbutton_mesh,'enable','on')


% --- Executes on button press in pushbutton_test.
function pushbutton_test_Callback(hObject, eventdata, handles)
global testV
testV=1;
run_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_prefs.
function pushbutton_prefs_Callback(hObject, eventdata, handles)
global f_canny_th1
global f_canny_th2
global f_canny_sig

global f_resize
global f_back
global f_histlim
global f_areamin
global f_areamax
global f_hmin
global f_hmin_split
global f_gstd
global f_pert_same
global f_frame_diff
global f_r_int

global f_relsz
global f_int_rej
global f_seg_min
global f_contmin
global f_acute
global f_pixelprox
global f_curve_std

global f_rect

global f_Nit
global f_kint

global v_method
global v_imtype

f_canny_th1 = str2num(get(handles.edit_cannylow,'String'));
f_canny_th2 = str2num(get(handles.edit_cannyhigh,'String'));
f_canny_sig = str2num(get(handles.edit_cannysigma,'String'));

f_resize = str2num(get(handles.edit_resize,'String'));
f_back = str2num(get(handles.edit_back,'String'));
f_histlim = str2num(get(handles.edit_histlim,'String'));
f_areamin = str2num(get(handles.edit_min_area,'String'));
f_areamax = str2num(get(handles.edit_max_area,'String'));
f_hmin = str2num(get(handles.edit_primary,'String'));
f_hmin_split = str2num(get(handles.edit_secondary,'String'));
f_gstd = str2num(get(handles.edit_smooth,'String'));
f_r_int = str2num(get(handles.edit_force_smooth,'String'));
f_pert_same = str2num(get(handles.edit_overlap,'String'));
f_frame_diff = str2num(get(handles.edit_frame,'String'));

v_indpt = get(handles.checkbox_indpt,'Value');
v_prox = get(handles.checkbox_prox,'Value');
v_save = get(handles.checkbox_save,'Value');
v_falsepos = get(handles.checkbox_falsepos,'Value');
v_advanced = get(handles.checkbox_advanced,'Value');
v_seed = get(handles.checkbox_seed,'Value');
v_manualtrack = get(handles.checkbox_manualtrack,'Value');
v_persist = get(handles.checkbox_mip_persist,'Value');
v_simplethres = get(handles.checkbox_simplethres,'Value');
v_colorcode = get(handles.checkbox_colorcode,'Value');
v_exclude = get(handles.checkbox_excludeedge,'Value');
v_nocontours = get(handles.checkbox_nocontours,'Value');


v_meshmotif = get(handles.checkbox_meshmotif,'Value');
v_mm_mesh = get(handles.checkbox_mm_mesh,'Value');
v_mt_mesh = get(handles.checkbox_mt_mesh,'Value');

[ftemp,ptemp] = uiputfile('Morphometrics_prefs.mat','Save parameters as ...');

if ftemp~=0
    save(fullfile(ptemp,ftemp),...
        'f_canny_th1',...
        'f_canny_th2',...
        'f_canny_sig',...
        'f_resize',...
        'f_back',...
        'f_histlim',...
        'f_areamin',...
        'f_areamax',...
        'f_hmin',...
        'f_hmin_split',...
        'f_gstd',...
        'f_r_int',...
        'f_pert_same',...
        'f_frame_diff',...
        'f_Nit',...
        'f_kint',...
        'f_relsz',...
        'f_int_rej',...
        'f_seg_min',...
        'f_contmin',...
        'f_acute',...
        'f_pixelprox',...
        'f_curve_std',...
        'f_rect',...
        'v_method',...
        'v_imtype',...
        'v_indpt',...
        'v_prox',...
        'v_save',...
        'v_falsepos',...
        'v_advanced',...
        'v_seed',...
        'v_manualtrack',...
        'v_persist',...
        'v_simplethres',...
        'v_colorcode',...
        'v_exclude',...
        'v_nocontours',...
        'v_meshmotif',...
        'v_mm_mesh',...
        'v_mt_mesh');
    set(handles.text_process,'String',['Parameter values saved to: ' ftemp])
end
set(handles.text_process,'String','Parameters saved.')
guidata(hObject,handles);


% --- Executes on button press in pushbutton_load_prefs.
function pushbutton_load_prefs_Callback(hObject, eventdata, handles)
global f_canny_th1
global f_canny_th2
global f_canny_sig

global f_resize
global f_back
global f_histlim
global f_areamin
global f_areamax
global f_hmin
global f_hmin_split
global f_gstd
global f_pert_same
global f_frame_diff
global f_r_int

global f_relsz
global f_int_rej
global f_seg_min
global f_contmin
global f_acute
global f_pixelprox
global f_curve_std

global f_rect

global f_Nit
global f_kint

global v_method
global v_imtype

[ftemp, ptemp] = uigetfile('*.mat', 'Choose Parameter File');

if ftemp~=0
    %load and update parameters
    try
        load(fullfile(ptemp,ftemp),'f_*','v_*')
        
        try
            set(handles.edit_resize,'String',num2str(f_resize));
        catch
            warning('Parameter ''f_resize'' missing in parameter file.')
        end
        try
            set(handles.edit_back,'String',num2str(f_back));
        catch
            warning('Parameter ''f_back'' missing in parameter file.')
        end
        try
            set(handles.edit_histlim,'String',num2str(f_histlim));
        catch
            warning('Parameter ''f_histlim'' missing in parameter file.')
        end
        try
            set(handles.edit_min_area,'String',num2str(f_areamin));
        catch
            warning('Parameter ''f_areamin'' missing in parameter file.')
        end
        try
            set(handles.edit_max_area,'String',num2str(f_areamax));
        catch
            warning('Parameter ''f_areamax'' missing in parameter file.')
        end
        try
            set(handles.edit_smooth,'String',num2str(f_gstd));
        catch
            warning('Parameter ''f_gstd'' missing in parameter file.')
        end
        try
            set(handles.edit_force_smooth,'String',num2str(f_r_int));
        catch
            warning('Parameter ''f_r_int'' missing in parameter file.')
        end
        try
            set(handles.edit_overlap,'String',num2str(f_pert_same));
        catch
            warning('Parameter ''f_pert_same'' missing in parameter file.')
        end
        try
            set(handles.edit_frame,'String',num2str(f_frame_diff));
        catch
            warning('Parameter ''f_frame_diff'' missing in parameter file.')
        end
        if ~exist('f_relsz','var')
            warning('Parameter ''f_relsz'' missing in parameter file.')
        end
        if ~exist('f_int_rej','var')
            warning('Parameter ''f_int_rej'' missing in parameter file.')
        end
        if ~exist('f_seg_min','var')
            warning('Parameter ''f_seg_min'' missing in parameter file.')
        end
        if ~exist('f_contmin','var')
            warning('Parameter ''f_contmin'' missing in parameter file.')
        end
        if ~exist('f_acute','var')
            warning('Parameter ''f_acute'' missing in parameter file.')
        end
        if ~exist('f_pixelprox','var')
            warning('Parameter ''f_pixelprox'' missing in parameter file.')
        end
        if ~exist('f_curve_std','var')
            warning('Parameter ''f_curve_std'' missing in parameter file.')
        end
        if ~exist('f_canny_th1','var')
            warning('Parameter ''f_canny_th1'' missing in parameter file.')
        end
        if ~exist('f_canny_th2','var')
            warning('Parameter ''f_canny_th2'' missing in parameter file.')
        end
        if ~exist('f_canny_sig','var')
            warning('Parameter ''f_canny_sig'' missing in parameter file.')
        end
        if ~exist('f_Nit','var')
            warning('Parameter ''f_Nit'' missing in parameter file.')
        end
        if ~exist('f_kint','var')
            warning('Parameter ''f_kint'' missing in parameter file.')
        end
        if ~exist('f_rect','var')
            warning('Parameter ''f_rect'' missing in parameter file.')
        end
        
        try
            switch v_method
                case 1
                    checkbox_seg_grad_Callback(hObject, eventdata, handles)
                case 2
                    checkbox_seg_lap_Callback(hObject, eventdata, handles)
                case 3
                    checkbox_seg_thres_Callback(hObject, eventdata, handles)
                case 4
                    checkbox_seg_canny_Callback(hObject, eventdata, handles)
            end
        catch
            warning('Parameters ''v_method'' missing in parameter file.')
        end
        
        try
            switch v_imtype
                case 1
                    checkbox_phase_Callback(hObject, eventdata, handles)
                case 2
                    checkbox_fluor_int_Callback(hObject, eventdata, handles)
                case 3
                    checkbox_fluor_periphery_Callback(hObject, eventdata, handles)
            end
        catch
            warning('Parameters ''v_imtype'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_indpt,'Value',v_indpt);
        catch
            warning('Parameters ''v_indpt'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_prox,'Value',v_prox);
        catch
            warning('Parameters ''v_prox'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_save,'Value',v_save);
        catch
            warning('Parameters ''v_save'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_falsepos,'Value',v_falsepos);
        catch
            warning('Parameters ''v_falsepos'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_advanced,'Value',v_advanced);
        catch
            warning('Parameters ''v_advanced'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_seed,'Value',v_seed);
        catch
            warning('Parameters ''v_seed'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_manualtrack,'Value',v_manualtrack);
        catch
            warning('Parameters ''v_manualtrack'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_mip_persist,'Value',v_persist);
        catch
            warning('Parameters ''v_persist'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_simplethres,'Value',v_simplethres);
        catch
            warning('Parameters ''v_simplethres'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_colorcode,'Value',v_colorcode);
        catch
            warning('Parameters ''v_colorcode'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_excludeedge,'Value',v_exclude);
        catch
            warning('Parameters ''v_exclude'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_nocontours,'Value',v_nocontours);
        catch
            warning('Parameters ''v_nccontours'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_meshmotif,'Value',v_meshmotif);
        catch
            warning('Parameter ''v_meshmotif'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_mm_mesh,'Value',v_mm_mesh);
        catch
            warning('Parameter ''v_mm_mesh'' missing in parameter file.')
        end
        
        try
            set(handles.checkbox_mt_mesh,'Value',v_mt_mesh);
        catch
            warning('Parameter ''v_mt_mesh'' missing in parameter file.')
        end
        
        load(fullfile(ptemp,ftemp),'f_hmin','f_hmin_split')
        try
            set(handles.edit_primary,'String',num2str(f_hmin));
        catch
            warning('Parameter ''f_hmin'' missing in parameter file.')
        end
        
        try
            set(handles.edit_secondary,'String',num2str(f_hmin_split));
        catch
            warning('Parameter ''f_hmin_split'' missing in parameter file.')
        end
    catch
        warning('The selected file either did not contain parameter values or is improperly formatted.')
    end
end
set(handles.text_process,'String','Parameters loaded.')
checkbox_advanced_Callback(hObject,eventdata,handles)
checkbox_prox_Callback(hObject,eventdata,handles)
guidata(hObject,handles);


% --- Executes on button press in pushbutton_mesh.
function pushbutton_mesh_Callback(hObject, eventdata, handles)
global mesh_motif_text
if get(handles.checkbox_meshmotif,'Value')
    %Set mesh file name processing motif
    prompt={'Enter the ''*CONTOURS*.mat'' batch file processing motif: '};
    name='Mesh Batch File Motif';
    numlines=1;
    defaultanswer={mesh_motif_text};
    
    mesh_motif_text_temp=char(inputdlg(prompt,name,numlines,defaultanswer));
    
    if isempty(mesh_motif_text_temp)
        return
    else
        mesh_motif_text=mesh_motif_text_temp;
    end
    
    %load mat data file to get path
    [file1,path1]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
    if file1==0
        return
    end
    
    %find files that match search string motif in current directory
    cd(path1)
    dir1=dir(mesh_motif_text);
    
    if isempty(dir1)
        warndlg('No files match the motif search string in this directory.', 'No files Match');
        return
    else
        for i=1:length(dir1)
            curr_file=dir1(i).name;
            curr_text=['Meshing ' curr_file ', ' num2str(i) ' of ' num2str(length(dir1))];
            set(handles.text_process,'String',curr_text)
            
            %select type of meshing
            if get(handles.checkbox_mm_mesh,'Value')     
                morphometrics_mesh(handles,curr_file,path1)
            else
                MT_mesh(handles,curr_file,path1)
            end
        end
    end
else
    %select type of meshing
    if get(handles.checkbox_mm_mesh,'Value')
        morphometrics_mesh(handles,[],[])
    else
        MT_mesh(handles,[],[])
    end
end


% --------------------------------------------------------------------
function adv_opts_Callback(hObject, eventdata, handles)
global f_relsz
global f_int_rej
global f_seg_min

global f_contmin
global f_acute
global f_pixelprox
global f_curve_std

prompt={'Diameter of relative noise filter (pxls)',...
    'False positive rejection parameter (>0)',...
    'Intensity histogram tail parameter (>0 and <0.45)',...
    'Minimum contour length (points)',...
    'Contour acute angle parameter(>0)',...
    'Cell proximity cutoff (pixels)',...
    'Standard deviation of curvature smoothing (pixels)',...
    'Lucky lotto numbers (pick 6)'};
name='Advanced Options';
numlines=1;
defaultanswer={num2str(f_relsz),num2str(f_int_rej),num2str(f_seg_min),num2str(f_contmin),num2str(f_acute),num2str(f_pixelprox),num2str(f_curve_std),...
    [num2str(round(100*rand)) ' ' num2str(round(100*rand)) ' ' num2str(round(100*rand)) ' ' num2str(round(100*rand)) ' ' num2str(round(100*rand)) ' ' num2str(round(100*rand))]};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer1=inputdlg(prompt,name,numlines,defaultanswer,options);

if ~isempty(answer1)
    temp1=str2num(answer1{1});
    if temp1>1
        f_relsz=temp1;
    else
        warndlg('The relative noise filter diameter must be greater than 1.','Invalid Option')
    end
    
    temp2=str2num(answer1{2});
    if temp2>0
        f_int_rej=temp2;
    else
        warndlg('The false positive rejection parameter must be greater than 0.','Invalid Option')
    end
    
    %{
    temp3=str2num(answer1{3});
    if temp3>=1
        f_seg_min=temp3;
    else
        warndlg('The minimum segmentation variance parmeter must be greater than or equal to 1.','Invalid Option')
    end
    %}
    
    temp3=str2num(answer1{3});
    if and(temp3>0,temp3<=0.45)
        f_seg_min=temp3;
    elseif and(temp3>0.45,temp3<=0.50)
        f_seg_min=temp3;
        warndlg('The intensity histogram tail weight parameter should be between 0 and 0.45, but the current value has been accepted.','Invalid Option')
    else
        warndlg('The intensity histogram tail weight parameter should be between 0 and 0.45.','Invalid Option')
    end
    
    temp4=str2num(answer1{4});
    if temp4>4
        f_contmin=temp4;
    else
        warndlg('The minimum number of contour points must be greater than 4.','Invalid Option')
    end
    
    temp5=str2num(answer1{5});
    if temp5>0
        f_acute=temp5;
    else
        warndlg('The acute angle removal parmeter must be greater than 0.','Invalid Option')
    end
    
    temp6=str2num(answer1{6});
    if temp6>1
        f_pixelprox=temp6;
    else
        warndlg('The cell proximity parmeter must be greater than 1 pixel.','Invalid Option')
    end
    
    temp7=str2num(answer1{7});
    if temp7>0
        f_curve_std=temp7;
    else
        warndlg('The curvature smoothing parmeter must be greater than 1 pixel.','Invalid Option')
    end
end


function checkbox_prox_Callback(hObject, eventdata, handles)
if get(handles.checkbox_prox,'Value')
    set(handles.edit_secondary,'enable','on')
else
    set(handles.edit_secondary,'enable','off')
end
guidata(hObject,handles);

function checkbox_indpt_Callback(hObject, eventdata, handles)
if get(handles.checkbox_indpt,'Value')
    set(handles.edit_overlap,'enable','off')
    set(handles.edit_frame,'enable','off')
    set(handles.checkbox_manualtrack,'enable','off')
    set(handles.checkbox_seed,'enable','off')
    set(handles.checkbox_seed,'Value',0)
else
    if get(handles.checkbox_advanced,'Value')
    set(handles.edit_overlap,'enable','on')
    set(handles.edit_frame,'enable','on')
    set(handles.checkbox_manualtrack,'enable','on')
    set(handles.checkbox_seed,'enable','on')
    checkbox_manualtrack_Callback(hObject, eventdata, handles)
    end
end
guidata(hObject,handles);

function checkbox_phase_Callback(hObject, eventdata, handles)
global v_imtype
v_imtype=1;
set(handles.checkbox_fluor_int,'Value',0)
set(handles.checkbox_fluor_periphery,'Value',0)
set(handles.checkbox_phase,'Value',1)
guidata(hObject,handles);

function checkbox_fluor_int_Callback(hObject, eventdata, handles)
global v_imtype
v_imtype=2;
set(handles.checkbox_phase,'Value',0)
set(handles.checkbox_fluor_periphery,'Value',0)
set(handles.checkbox_fluor_int,'Value',1)
guidata(hObject,handles);

function checkbox_fluor_periphery_Callback(hObject, eventdata, handles)
global v_imtype
global v_method
v_imtype=3;

if v_method==2
    checkbox_phase_Callback(hObject, eventdata, handles)
    warndlg('Laplacian segmentation is incompatible with Fluorescence (peripheral) images.', 'Method Mistmatch');
else
    set(handles.checkbox_phase,'Value',0)
    set(handles.checkbox_fluor_int,'Value',0)
    set(handles.checkbox_fluor_periphery,'Value',1)
end

if v_method==3
    checkbox_fluor_int_Callback(hObject, eventdata, handles)
    warndlg('Adaptive segmentation is incompatible with Fluorescence (peripheral) images.', 'Method Mistmatch');
else
    set(handles.checkbox_phase,'Value',0)
    set(handles.checkbox_fluor_int,'Value',0)
    set(handles.checkbox_fluor_periphery,'Value',1)
end
guidata(hObject,handles);

% --- Executes on button press in checkbox_seed.
function checkbox_seed_Callback(hObject, eventdata, handles)
if get(handles.checkbox_seed,'Value')
    set(handles.checkbox_manualtrack,'enable','off')
    set(handles.checkbox_manualtrack,'Value',0)
    set(handles.edit_overlap,'enable','off')
    set(handles.edit_frame,'enable','off')
    set(handles.text_process,'String','Images should be aligned (Options menu) when seeding contours.')
else
    set(handles.checkbox_manualtrack,'enable','on')
    set(handles.edit_overlap,'enable','on')
    set(handles.edit_frame,'enable','on')
end

function checkbox_help_Callback(hObject, eventdata, handles)
if get(handles.checkbox_fluor_periphery,'Value')
    warndlg('This feature is not compatible with peripheral fluorescence.', 'Feature Unavailable');
    set(handles.checkbox_help,'Value',0)
end

function edit_min_area_Callback(hObject, eventdata, handles)
global f_areamin
global f_areamax

f_areamin=round(str2num(get(handles.edit_min_area,'String')));
if f_areamin<9
    f_areamin=9;
end
if f_areamin>f_areamax
    f_areamax=2*f_areamin;
end
set(handles.edit_min_area,'String',num2str(f_areamin))
set(handles.edit_max_area,'String',num2str(f_areamax))
guidata(hObject,handles);

function edit_max_area_Callback(hObject, eventdata, handles)
global f_areamin
global f_areamax

if or(strcmp(get(handles.edit_max_area,'String'),'Inf'),strcmp(get(handles.edit_max_area,'String'),'inf'))
    set(handles.edit_max_area,'String','1000000')
    errordlg('The maximum object size must be finite.', 'Max size error');
end

f_areamin=round(str2num(get(handles.edit_min_area,'String')));
f_areamax=round(str2num(get(handles.edit_max_area,'String')));

if f_areamax<=f_areamin
    f_areamax=f_areamin+10;
end

set(handles.edit_max_area,'String',num2str(f_areamax))
guidata(hObject,handles);

function edit_smooth_Callback(hObject, eventdata, handles)
global f_gstd

f_gstd=str2num(get(handles.edit_smooth,'String'));

if f_gstd<0
    f_gstd=0;
end
set(handles.edit_smooth,'String',num2str(f_gstd))
guidata(hObject,handles);

function edit_primary_Callback(hObject, eventdata, handles)
global f_hmin
global v_method

f_hmin=str2num(get(handles.edit_primary,'String'));

if and(or(v_method==2,v_method==1),or(f_hmin<0,f_hmin>1))
    f_hmin=0.1;
end

if all([v_method==3,get(handles.checkbox_simplethres,'Value'),or(f_hmin<0,f_hmin>1)])
    f_hmin=0.1;
end

set(handles.edit_primary,'String',num2str(f_hmin))
guidata(hObject,handles);

function edit_secondary_Callback(hObject, eventdata, handles)
global f_hmin_split

f_hmin_split=str2num(get(handles.edit_secondary,'String'));

if f_hmin_split<1
    f_hmin_split=1;
end

set(handles.edit_secondary,'String',num2str(f_hmin_split))
guidata(hObject,handles);

function edit_test_frame_Callback(hObject, eventdata, handles)
global path1
global fnames

if path1==0
    warndlg('First load an input file, then adjust this parameter.', 'Load Data First');
    return
end

%get image info
Im_info=imfinfo(fullfile(path1,fnames(1).name));

t_frame=round(str2num(get(handles.edit_test_frame,'String')));

if or(t_frame<1,t_frame>length(Im_info))
    warndlg('This frame number is out of range for this stack.', 'Frame Out of Range');
    set(handles.edit_test_frame,'String',num2str(1))
    return
end
guidata(hObject,handles);

function edit_force_smooth_Callback(hObject, eventdata, handles)
global f_r_int

f_r_int=str2num(get(handles.edit_force_smooth,'String'));

if f_r_int<=0
    f_r_int=0;
end

set(handles.edit_force_smooth,'String',num2str(f_r_int))
guidata(hObject,handles);

function edit_overlap_Callback(hObject, eventdata, handles)
global f_pert_same

f_pert_same=str2num(get(handles.edit_overlap,'String'));

if or(f_pert_same<=0,f_pert_same>=1)
    f_pert_same=0.6;
end

if f_pert_same<0.5
    helpdlg('Setting the overlap below 0.5 in dense cell populations may cause loss of cell data.','Just so you know.');
end

set(handles.edit_overlap,'String',num2str(f_pert_same))
guidata(hObject,handles);

function edit_frame_Callback(hObject, eventdata, handles)
global f_frame_diff
global path1
global fnames

if path1==0
    warndlg('First load an input file, then change this value.', 'Load Data First');
    set(handles.edit_frame,'String',num2str(f_frame_diff))
    return
end

%get image info
Im_info=imfinfo(fullfile(path1,fnames(1).name));

f_frame_diff=str2num(get(handles.edit_frame,'String'));

if f_frame_diff<1
    f_frame_diff=1;
elseif f_frame_diff>length(Im_info)
    f_frame_diff=length(Im_info);
end

set(handles.edit_frame,'String',num2str(f_frame_diff))
guidata(hObject,handles);

function edit_histlim_Callback(hObject, eventdata, handles)
global f_histlim

f_histlim=str2num(get(handles.edit_histlim,'String'));

if or(f_histlim<0,f_histlim>=0.5)
    f_histlim=0;
    set(handles.edit_histlim,'String',num2str(f_histlim))
    warndlg('The historgram limit must be greater than or equal to 0 and less than 0.5, and probably close to zero.','Historgram Limit Error')
end

guidata(hObject,handles);

function edit_back_Callback(hObject, eventdata, handles)
global f_back

f_back=str2num(get(handles.edit_back,'String'));
if f_back<0
    f_back=0;
    set(handles.edit_back,'String',num2str(f_back))
end
guidata(hObject,handles);

function edit_resize_Callback(hObject, eventdata, handles)
global f_resize

f_resize=str2num(get(handles.edit_resize,'String'));
if or(f_resize<0.1,f_resize>10)
    f_resize=1;
    set(handles.edit_resize,'String','1')
end

function checkbox_mt_mesh_Callback(hObject, eventdata, handles)
set(handles.checkbox_mm_mesh,'Value',0)
set(handles.checkbox_mt_mesh,'Value',1)

function checkbox_mm_mesh_Callback(hObject, eventdata, handles)
set(handles.checkbox_mm_mesh,'Value',1)
set(handles.checkbox_mt_mesh,'Value',0)

function checkbox_manualtrack_Callback(hObject, eventdata, handles)
if get(handles.checkbox_manualtrack,'Value')
    set(handles.edit_overlap,'enable','off')
    set(handles.edit_frame,'enable','off')
    set(handles.checkbox_seed,'enable','off')
    set(handles.checkbox_seed,'Value',0)
else
    set(handles.edit_overlap,'enable','on')
    set(handles.edit_frame,'enable','on')
    set(handles.checkbox_seed,'enable','on')
end

% --------------------------------------------------------------------
function checkbox_seg_grad_Callback(hObject, eventdata, handles)
global f_int_rej
f_int_rej=1;

global v_method
v_method=1;

set(handles.checkbox_seg_grad,'Value',1)
set(handles.checkbox_seg_lap,'Value',0)
set(handles.checkbox_seg_thres,'Value',0)
set(handles.checkbox_seg_canny,'Value',0)
set(handles.checkbox_mip_persist,'enable','off')
set(handles.checkbox_simplethres,'enable','off')
set(handles.edit_cannylow,'enable','off')
set(handles.edit_cannyhigh,'enable','off')
set(handles.edit_cannysigma,'enable','off')
set(handles.text42,'enable','off')
set(handles.text43,'enable','off')
set(handles.text44,'enable','off')

set(handles.edit_histlim,'enable','on')

set(handles.checkbox_falsepos,'Value',1)
set(handles.checkbox_falsepos,'enable','on')

set(handles.checkbox_prox,'Value',0)
set(handles.checkbox_prox,'enable','on')

set(handles.edit_secondary,'enable','off')
set(handles.edit_primary,'enable','on')
set(handles.text23,'enable','on')

set(handles.checkbox_help,'enable','on')

set(handles.text23,'String','Primary threshold (hmin)')
set(handles.text23,'TooltipString','Watershed thresholding parameter for initial segmentation.');

checkbox_advanced_Callback(hObject, eventdata, handles)
edit_primary_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function checkbox_seg_lap_Callback(hObject, eventdata, handles)
global f_int_rej
f_int_rej=3;

global v_method
v_method=2;

set(handles.checkbox_seg_grad,'Value',0)
set(handles.checkbox_seg_lap,'Value',1)
set(handles.checkbox_seg_thres,'Value',0)
set(handles.checkbox_seg_canny,'Value',0)
set(handles.checkbox_mip_persist,'enable','off')
set(handles.checkbox_simplethres,'enable','off')
set(handles.edit_cannylow,'enable','off')
set(handles.edit_cannyhigh,'enable','off')
set(handles.edit_cannysigma,'enable','off')
set(handles.text42,'enable','off')
set(handles.text43,'enable','off')
set(handles.text44,'enable','off')

set(handles.edit_histlim,'enable','off')

set(handles.checkbox_falsepos,'Value',1)
set(handles.checkbox_falsepos,'enable','on')

set(handles.checkbox_prox,'Value',0)
set(handles.checkbox_prox,'enable','on')

set(handles.edit_secondary,'enable','on')
set(handles.edit_primary,'enable','on')
set(handles.text23,'enable','on')

set(handles.checkbox_help,'enable','off')

set(handles.text23,'String','Primary threshold (hmin)')
set(handles.text23,'TooltipString','Watershed thresholding parameter for initial segmentation.');

if get(handles.checkbox_fluor_periphery,'Value')
    warndlg('Laplacian Segmentation does not work with peripheral fluorescence.', 'Method mismatch');
    checkbox_seg_grad_Callback(hObject, eventdata, handles)
end

checkbox_advanced_Callback(hObject, eventdata, handles)
edit_primary_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function checkbox_seg_thres_Callback(hObject, eventdata, handles)
global v_method
v_method=3;

set(handles.checkbox_seg_grad,'Value',0)
set(handles.checkbox_seg_lap,'Value',0)
set(handles.checkbox_seg_thres,'Value',1)
set(handles.checkbox_seg_canny,'Value',0)
set(handles.checkbox_mip_persist,'enable','on')
set(handles.checkbox_simplethres,'enable','on')
set(handles.edit_cannylow,'enable','off')
set(handles.edit_cannyhigh,'enable','off')
set(handles.edit_cannysigma,'enable','off')
set(handles.text42,'enable','off')
set(handles.text43,'enable','off')
set(handles.text44,'enable','off')

set(handles.edit_histlim,'enable','on')

set(handles.checkbox_falsepos,'Value',1)
set(handles.checkbox_falsepos,'enable','on')

set(handles.checkbox_prox,'Value',0)
set(handles.checkbox_prox,'enable','on')

set(handles.edit_secondary,'enable','off')
set(handles.edit_primary,'enable','on')
set(handles.text23,'enable','on')
if get(handles.checkbox_simplethres,'Value')
    set(handles.text23,'String','Intensity threshold')
else
    set(handles.text23,'String','Intensity threshold (STD)')
end

set(handles.checkbox_help,'enable','on')

set(handles.text23,'TooltipString','The number of standard deviations above the selected background standard deviation that triggers segmentation of a pixel region.');

set(handles.edit_primary,'String','3')

checkbox_advanced_Callback(hObject, eventdata, handles)
edit_primary_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function checkbox_seg_canny_Callback(hObject, eventdata, handles)
global f_canny_th1
global f_canny_th2
global f_canny_sig
global v_method
v_method=4;

set(handles.checkbox_seg_grad,'Value',0)
set(handles.checkbox_seg_lap,'Value',0)
set(handles.checkbox_seg_thres,'Value',0)
set(handles.checkbox_seg_canny,'Value',1)
set(handles.checkbox_mip_persist,'enable','off')
set(handles.checkbox_simplethres,'enable','off')
set(handles.edit_cannylow,'enable','on')
set(handles.edit_cannyhigh,'enable','on')
set(handles.edit_cannysigma,'enable','on')
set(handles.text42,'enable','on')
set(handles.text43,'enable','on')
set(handles.text44,'enable','on')

set(handles.edit_histlim,'enable','off')

set(handles.checkbox_falsepos,'Value',0)
set(handles.checkbox_falsepos,'enable','off')

set(handles.checkbox_prox,'Value',0)
set(handles.checkbox_prox,'enable','on')

set(handles.edit_secondary,'enable','on')
set(handles.edit_primary,'enable','off')
set(handles.text23,'enable','off')

set(handles.checkbox_help,'enable','off')

set(handles.text23,'String','Primary threshold (hmin)');
set(handles.text23,'TooltipString','Watershed thresholding parameter for initial segmentation.');

if get(handles.checkbox_fluor_periphery,'Value')
    warndlg('Canny Segmentation does not work with peripheral fluorescence.', 'Method mismatch');
    checkbox_seg_grad_Callback(hObject, eventdata, handles)
end

checkbox_advanced_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox_advanced.
function checkbox_advanced_Callback(hObject, eventdata, handles)
if get(handles.checkbox_advanced,'Value')
    set(handles.edit_resize,'enable','on')
    set(handles.text_resize,'enable','on')
    set(handles.edit_back,'enable','on')
    set(handles.text_back,'enable','on')
    set(handles.edit_histlim,'enable','on')
    set(handles.text_histlim,'enable','on')
    set(handles.edit_smooth,'enable','on')
    set(handles.text_smooth,'enable','on')
    set(handles.edit_force_smooth,'enable','on')
    set(handles.text_force_smooth,'enable','on')
    set(handles.edit_overlap,'enable','on')
    set(handles.text_overlap,'enable','on')
    set(handles.edit_frame,'enable','on')
    set(handles.text_frame,'enable','on')
    if ~get(handles.checkbox_indpt,'Value')
        set(handles.checkbox_manualtrack,'enable','on')
        set(handles.checkbox_seed,'enable','on')
    else
        set(handles.checkbox_manualtrack,'enable','off')
        set(handles.checkbox_seed,'enable','off')
    end
    set(handles.checkbox_seg_canny,'enable','on')
    if get(handles.checkbox_seg_canny,'Value')
        set(handles.edit_cannylow,'enable','on')
        set(handles.edit_cannyhigh,'enable','on')
        set(handles.edit_cannysigma,'enable','on')
        set(handles.text42,'enable','on')
        set(handles.text43,'enable','on')
        set(handles.text44,'enable','on')
    else
        set(handles.edit_cannylow,'enable','off')
        set(handles.edit_cannyhigh,'enable','off')
        set(handles.edit_cannysigma,'enable','off')
        set(handles.text42,'enable','off')
        set(handles.text43,'enable','off')
        set(handles.text44,'enable','off')
    end
else
    set(handles.edit_resize,'enable','off')
    set(handles.text_resize,'enable','off')
    set(handles.edit_back,'enable','off')
    set(handles.text_back,'enable','off')
    set(handles.edit_histlim,'enable','off')
    set(handles.text_histlim,'enable','off')
    set(handles.edit_smooth,'enable','off')
    set(handles.text_smooth,'enable','off')
    set(handles.edit_force_smooth,'enable','off')
    set(handles.text_force_smooth,'enable','off')
    set(handles.edit_overlap,'enable','off')
    set(handles.text_overlap,'enable','off')
    set(handles.edit_frame,'enable','off')
    set(handles.text_frame,'enable','off')
    set(handles.checkbox_manualtrack,'enable','off')
    set(handles.checkbox_seed,'enable','off')
    set(handles.checkbox_seg_canny,'enable','off')
    set(handles.edit_cannylow,'enable','off')
    set(handles.edit_cannyhigh,'enable','off')
    set(handles.edit_cannysigma,'enable','off')
    set(handles.text42,'enable','off')
    set(handles.text43,'enable','off')
    set(handles.text44,'enable','off')
    
    set(handles.checkbox_manualtrack,'Value',0)
    set(handles.checkbox_seed,'Value',0)
    if get(handles.checkbox_seg_canny,'Value')
        set(handles.checkbox_seg_grad,'Value',1)
    end
    set(handles.checkbox_seg_canny,'Value',0)    
end

% --------------------------------------------------------------------
function menu_imagealign_Callback(hObject, eventdata, handles)
stack_align

% --------------------------------------------------------------------
function menu_combine_Callback(hObject, eventdata, handles)
stack_maker

% --------------------------------------------------------------------
function menu_create_cell_table_Callback(hObject, eventdata, handles)
invert_frame(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_view_contours_Callback(hObject, eventdata, handles)
view_contours

% --------------------------------------------------------------------
% function menu_update_Callback(hObject, eventdata, handles) %Updates will
% be through github not the software directly

% %version control check
% try
%     persist_url='https://www.dropbox.com/s/ul5djjt8krl2kei/morphometrics_version.txt?raw=1';
%     ver_num='v1.1';  %<--- update latest version number here
%     display('Attempting connection with update server at:')
%     display([' ' persist_url])
%     vers_test=urlread(persist_url,'Timeout',3);
%     try
%         if ~strcmp(ver_num,vers_test)
%             v0=questdlg('A new version of Morphometrics is available, would you like to download it?', ...
%                 'Update Morphometrics','Yes', 'No', 'Yes');
%             if strcmp(v0,'Yes')
%                 web('https://www.dropbox.com/s/7b48ppqs1pcyq3t/Morphometrics_v1.1.zip?dl=1')
%             end
%         else
%             helpdlg(['Your version of Morphometrics is up-to-date! (' ver_num ')'],'Up to date')
%         end
%     catch
%         warndlg('An error occurred while attemping to download the latest version of Morphometrics.','Download Error');
%     end
% catch
%     warndlg('Morphometrics could not connect to the update server.','Update Error')
% end

% --- Executes on button press in pushbutton_roi.
function pushbutton_roi_Callback(hObject, eventdata, handles)
global path1
global fnames
global h1
global X_roi
global Y_roi

if path1==0
    warndlg('First load an input file, then try drawing an ROI.', 'Load Data First');
    return
end

%delete previous handles
try temp2=getPosition(h1);
    X_roi=0;
    Y_roi=0;
    delete(h1)
    drawnow
catch
end

%get current test frame
h0 = imfreehand(handles.axes1,'Closed',1);
if ~isempty(h0)
    setClosed(h0,'True')
    
    %get coordinates
    temp1=getPosition(h0);
    X_roi=temp1(:,1);
    Y_roi=temp1(:,2);
    
    %fix out of bounds
    Im_info=imfinfo(fullfile(path1,fnames(1).name));
    
    X_roi(X_roi<1)=1;
    X_roi(X_roi>Im_info(1).Width)=Im_info(1).Width;
    Y_roi(Y_roi<1)=1;
    Y_roi(Y_roi>Im_info(1).Height)=Im_info(1).Height;
    
    h1=h0;
end

% --- Executes on button press in pushbutton_deleteroi.
function pushbutton_deleteroi_Callback(hObject, eventdata, handles)
global X_roi
global Y_roi
global h1

X_roi=0;
Y_roi=0;

try temp2=getPosition(h1);
    delete(h1)
    drawnow
catch
end

% --- Executes on button press in pushbutton_exit.
function pushbutton_exit_Callback(hObject, eventdata, handles)
set(handles.text_process,'String','Good bye!')
display('Good bye!')
drawnow()
clearvars
close('Morphometrics version 2.0')

% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
global stopV
global testV
stopV=1;
testV=0;

set(handles.text_process,'String','Processing stopped by user ...')
title('Processing stopped by user ...')

enable_buttons(handles)


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
try
    open('Morphometrics_v1.1_Help.pdf')
catch
    errordlg('Could not locate Morphometrics help file.')
end
%popupmessage('Morphometrics_help.txt','Morphometrics Help')


%*************************************************************************
%******    NON-USER FUNCTIONS ********************************************
%*************************************************************************

% --- Outputs from this function are returned to the command line.
function varargout = morphometrics_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function text_motif_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_save.
function checkbox_save_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_static_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text16_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function checkbox_prox_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function checkbox_save_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function open_mult_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_motif_static_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function checkbox_fluor_periphery_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function checkbox_phase_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function checkbox_indpt_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_force_smooth_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_motif_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_min_area_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_frame_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_overlap_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_test_frame_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_secondary_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_primary_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_smooth_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_max_area_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_back_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_histlim_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_falsepos.
function checkbox_falsepos_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_meshmotif.
function checkbox_meshmotif_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_resize_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_savetest.
function checkbox_savetest_Callback(hObject, eventdata, handles)

function menu_file_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_keeptemp.
function checkbox_keeptemp_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_excludeedge.
function checkbox_excludeedge_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_keeptemp.
function checkbox_colorcode_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_options_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_enable_buttons_Callback(hObject, eventdata, handles)
enable_buttons(handles)


% --------------------------------------------------------------------
function menu_addpath_Callback(hObject, eventdata, handles)
global home_path
addpath(home_path)
addpath(fullfile(home_path,'Morphometrics_CL'))
savepath;

%turn off path add question dialog
fid= fopen(fullfile(home_path,'ask_add_to_path.txt'),'w+');
fprintf(fid,'%u',0);
fclose(fid);

disp(['Added to path:  ' home_path]);
disp(['Added to path:  ' fullfile(home_path,'Morphometrics_CL')]);


% --------------------------------------------------------------------
% function menu_noupdate_Callback(hObject, eventdata, handles)
% global home_path
% %turn off version control
% fid= fopen(fullfile(home_path,'update_on_start.txt'),'w+');
% fprintf(fid,'%u',0);
% fclose(fid);


% --- Executes on button press in checkbox_mip_persist.
function checkbox_mip_persist_Callback(hObject, eventdata, handles)
set(handles.checkbox_simplethres,'Value',0)
if get(handles.checkbox_simplethres,'Value')
    set(handles.text23,'String','Intensity threshold')
else
    set(handles.text23,'String','Intensity threshold (STD)')
end

% --- Executes on button press in checkbox_simplethres.
function checkbox_simplethres_Callback(hObject, eventdata, handles)
set(handles.checkbox_mip_persist,'Value',0)
if get(handles.checkbox_simplethres,'Value')
    set(handles.text23,'String','Intensity threshold')
else
    set(handles.text23,'String','Intensity threshold (STD)')
end


function edit_cannylow_Callback(hObject, eventdata, handles)
global f_canny_th1
temp1=str2num(get(handles.edit_cannylow,'String'));

if temp1<0
    temp1=1;
elseif temp1>str2num(get(handles.edit_cannyhigh,'String'))
    temp1=str2num(get(handles.edit_cannylow,'String'))-0.1;
    if temp1<0
        temp1=0;
    end
end

f_canny_th1=temp1;
set(handles.edit_cannylow,'String',num2str(f_canny_th1))
% --- Executes during object creation, after setting all properties.
function edit_cannylow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cannyhigh_Callback(hObject, eventdata, handles)
global f_canny_th2
temp1=str2num(get(handles.edit_cannyhigh,'String'));

if temp1>1
    temp1=1;
elseif temp1<str2num(get(handles.edit_cannylow,'String'))
    temp1=str2num(get(handles.edit_cannylow,'String'))+0.1;
    if temp1>1
        temp1=1;
    end
end

f_canny_th2=temp1;
set(handles.edit_cannyhigh,'String',num2str(f_canny_th2))
% --- Executes during object creation, after setting all properties.
function edit_cannyhigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cannysigma_Callback(hObject, eventdata, handles)
global f_canny_sig
temp1=str2num(get(handles.edit_cannysigma,'String'));

if temp1<0
    temp1=0;
end

f_canny_sig=temp1;
set(handles.edit_cannysigma,'String',num2str(f_canny_sig))
% --- Executes during object creation, after setting all properties.
function edit_cannysigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_dataprep_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_analysis_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_shapestats_Callback(hObject, eventdata, handles)
intensimetrics
disp(['Launching intensity analysis scripts ...'])


% --------------------------------------------------------------------
function quality_check_Callback(hObject, eventdata, handles)
quality_check


% --- Executes on button press in checkbox_nocontours.
function checkbox_nocontours_Callback(hObject, eventdata, handles)

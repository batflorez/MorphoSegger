
function morphometrics_mask_cl(stackname,paramname,varargin)
% Tristan Ursell, 2014.
% Modified by Andres Florez 04/20/2020

% Example:
% morphometrics_cl('mask.tif','Morphometrics_prefs_mask_CL.mat',[],1) with
% meshing -  Andres Florez 04/20/20 This version has been modified to work
% on binary masks by deleting saturarion detection

% stackname = full path to the image stack being processed
% paramname = full path to the parameter file
% rect1 = coordinates of a box for background determination during adaptive
% threshold segmentation [X0,Y0,Wid,Hgt] (NOT USED)
% meshq = a binary variable to tell morphometrics to perform meshing
%



%parse varargin input
if ~isempty(varargin)
    if length(varargin)==1
        rect1=varargin{1};
        meshq=0;
    else
        rect1=varargin{1};
        meshq=varargin{2};
    end
else
    rect1=0;
    meshq=0;
end

%close any open files
fclose all;

%initialize all global variables
global f_resize %image resizing factor
global f_back %background removal filter size
global f_histlim %histogram percentage removal for contrasting
global f_areamin %minimum object area
global f_areamax %maximum object area
global f_hmin %primary watershed threshold
global f_hmin_split %secondary pixel distance threshold
global f_gstd %initial Gaussian smoothing filter STD
global f_r_int %Gaussian force smoothing STD
global f_pert_same %percent overlap for cells to be related
global f_frame_diff %number of frames in back/front to look for related cells

%segmentation specific
global f_canny_th1 %lower canny threshold, triggers edge start
global f_canny_th2 %upper canny threshold, triggers edge end
global f_canny_sig %canny sigma factor, determines smoothness

%advanced options
global f_relsz %relative noise filter disk size (diameter)
global f_int_rej %false positive rejection parameter
global f_seg_min %minimum relative variacne required for any segmentation
global f_contmin %minimum number of points a contour may have
global f_acute %parameter that controls removal of acute angles in contours
global f_pixelprox %number of pixels below which two cells are considered proximal
global f_curve_std %standard deviation (in pixels) of guassian convolution curvature filter
global v_simplethres %Simple threshold true -This parameter was added back- Andres Florez 05/06/20

%contour fitting variables
global f_Nit %number of contour fitting iterations
global f_kint %rate of contour fitting

%methods and image type
global v_method %sets the method variable (gradient, laplacian, adaptive threshold)
global v_imtype %sets type of image being used (phase dark, interior fluor, peripheral fluor)
global v_seed %select whether to use seeded contours
global v_falsepos %select whether to reject false positives
global v_prox %select whether cells are in proximity
global v_indpt %select whether images are independent
global v_mm_mesh %use morphometrics full mesh and branching
global v_mt_mesh %use pill 1D mesh

%load(paramname,'f_*','v_*') % Uncomment this if you want to run it without
% the pipeline
v2struct(paramname)  % Unpack morphometrics parameters from structure

%check for presence of all inputs
q=0;
paramlist={       
    'f_resize',...
    'f_back',...
    'f_histlim',...
    'f_areamin',...
    'f_areamax',...
    'f_hmin',...
    'f_hmin_split',...
    'f_gstd',...
    'f_r_int',...
    'f_Nit',...
    'f_kint',...
    'f_pert_same',...
    'f_frame_diff',...
    'f_relsz',...
    'f_int_rej',...
    'f_seg_min',...
    'f_contmin',...
    'f_acute',...
    'f_pixelprox',...
    'f_curve_std',...
    'f_canny_th1',...
    'f_canny_th2',...
    'f_canny_sig',...
    'v_method',...
    'v_imtype',...
    'v_indpt',...
    'v_simplethres',... %This was added back - Andres Florez 05/06/20
    'v_prox',...
    'v_falsepos',...
    'v_seed',...
    'v_mm_mesh',...
    'v_mt_mesh'};
for i=1:length(paramlist)
    if exist(paramlist{i},'var')
        if isempty(eval(paramlist{i}))
            disp(['Parameter file is missing ' '' paramlist{i} '''.'])
            q=1;
        end
    else
        disp(['Parameter file is missing ' '' paramlist{i} '''.'])
        q=1;
    end
end
if q
    disp('Morphometrics cannot run without these variables.')
    return
end

% %check method / image compatibility
% if and(v_method==3,v_imtype==3)
%     %does not work with flourescence periphery
%     disp('Adaptive Threshold segmentation is incompatible with periphery images.', ' ');
%     return
% end

tic
Im_info=imfinfo(stackname);
Nim=length(Im_info);

%get file name parts
[pathstr, name, ext] = fileparts(stackname);
basename=fullfile(pathstr,name);

%get rid of temporary files if they exist
dir1=dir([basename '_G*.tif']);
if ~isempty(dir1)
    delete([basename '_Gparent.tif']);
    delete([basename,'_Gsegt.tif']);
end

disp(['Starting segmentation of: ' strtrim(stackname)])

%% LOOP OVER INDIVIDUAL IMAGES

for j=1:Nim
    Im0=imread(stackname,j);
    
%     %perform saturation detection  - This is not necessary for binary masks -Andres Florez 04/25/20
%     if Im_info(1).BitDepth==8
%         sat_test=or(Im0==0,Im0==(2^8-1));
%     else
%         sat_test=or(Im0==0,Im0==(2^16-1));
%     end
%     
%     sat_props=regionprops(sat_test,'Area');
%     
%     if any([sat_props.Area]>=9)
%         disp(['Warning: saturation detected on frame ' num2str(j) '.'])
%     end
%     

%%     SEGMENTATION
 
    %segment the images based on binary masks from Supersegger
    
    
    mask1=ones(size(Im0));
    %In case inverting the image is necessary - Andres Florez 04/20/20
    %Im0=logical(Im0); 
    %Im0=not(Im0);
    [segdata,Imag_out,Iseg_out,~]=simply_segment_cl(Im0,mask1,j,rect1);
    
    %save output segmentation data to main data structure
    if ~isempty(segdata)
        frame(j)=segdata;
    elseif and(j==1,0)
        frame(j).num_objs=frame(1).num_objs;
    else
        frame(j).num_objs=0;
    end
    
    %save output to stack for later processing
    image_save(Iseg_out,[basename '_Gparent.tif']) 
    image_save(Imag_out,[basename '_Gsegt.tif']) 
    
    disp(['Frame ' num2str(j) ' of ' num2str(Nim) ', ' num2str(frame(j).num_objs) ' objects'])
    
    clear Imag_out 
    clear Iseg_out
    clear segdata
    %-------------------------------
    %-------------------------------
end

%check to make sure something was found
if sum([frame.num_objs])<1
    disp('No objects were found in this stack using current parameter values. Consider adjusting the parameter values.')
    return
end


%%     LINEAGE

if all([~v_indpt,~v_seed,Nim>1])
    %find unique cells in stack
    try
        cell=family_cl([basename '_Gparent.tif'],f_pert_same,0,f_frame_diff);
        
        %update frame data structure to contain cell lineage data
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
    catch er1
        disp('An error was encountered while trying to create the cell lineage.')
        disp('Try adjusting parameter values, or try independent images.')
        disp(' ')
        disp('Analysis will continue using independent cell analysis on this stack.')
        disp(' ')
        disp(er1)
        
        clear cell
    end
end

%update frame data structure to contain *seeded* cell lineage
if v_seed
    for j=1:Nim
        frame(j).num_objs=frame(1).num_objs;
        for k=1:frame(1).num_objs
            frame(j).object(k).bw_label=k;
            frame(j).object(k).cellID=k;
        end
    end
end

%give cells unique names when frames are independent
if or(v_indpt,Nim==1)
    q=0;
    for j=1:Nim
        for k=1:frame(j).num_objs
            q=q+1;
            
            frame(j).object(k).bw_label=k;
            frame(j).object(k).cellID=q;
        end
    end
end

%determine maximum cell number
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
%-------------------------------
%-------------------------------


%% PROXIMITY TESTING & CONTOUR FITTING

disp('Fitting contours ...')
for j=1:Nim
    if and(frame(j).num_objs~=0,isfield(frame(j).object,'cellID'))
        temp_frame=frame(j);
        if and(j>1,v_seed)
            temp_frame_m1=frame(j-1);
            [frameout]=fit_contours_cl(basename,temp_frame,j,Nim,temp_frame_m1);
        else
            [frameout]=fit_contours_cl(basename,temp_frame,j,Nim);
        end
        frame(j)=frameout;
    end
end

%% CREATE CELL ID TABLE

% function was extracted from intensimetrics.m -  Andres Florez

cells=struct;
for i=1:Ncell
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
cells_ID_table =cells;
%% Deleting temporary files

%get rid of temporary files if they exist
dir1=dir([basename '_G*.tif']);
if ~isempty(dir1)
    delete([basename '_Gparent.tif']);
    delete([basename '_Gsegt.tif']);
end

%% OUTPUT FILES (move this to the morphometrics folder)

%assemble full data structure and write output data
outname=[basename '_' date '_CONTOURS.mat'];
save(outname,'cells_ID_table','frame','outname','Ncell','v_*','f_*')
disp('All relevant data has been saved to the file:')
disp([outname])

%% MESHING

if meshq
    if v_mm_mesh
        disp('Starting Voronoi meshing ...')
        morphometrics_mesh_cl(outname)
        disp('Meshing complete.')
    end
    if v_mt_mesh
        disp('Starting 1D meshing ...')
        MT_mesh_cl(outname)
        disp('Meshing complete.')
    end
end

t1=toc;
disp(['Finished in ' num2str(round(10*t1/60)/10) ' minutes.'])
clearvars


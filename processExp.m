%%%%%%%%%%%%% processExp Script - Pipeline for Imaging Analysis %%%%%%%%%%%
%                                                                         %
%     @@@@@@@@@@@@@@@@@@          @   @        @@@@@@@@@@@@@@@@@@@@       %
%         @@@@@@@@@@@@@@@@@       @@@@@       @@@@@@@@@@@@@@@@@           %
%             @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@              %
%               @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                %
%                 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                 %
%                          @@@@@@@@@@@@@@@@@@@                            %
%                               @@@@@@@@@                                 %
%                                  @@@@                                   %
%                                   @                                     %
%                                                                         %
% Process ND2 files, split channels, converts to tif, performs segmentation
% with SuperSegger and analysis with Morphometrics v2. . The current
% script was modified from SuperSegger. Copyright (C) 2016 Wiggins Lab.

% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

%%          (MODIFY THIS FILE AND SAVE IT IN DATA FOLDER)                %%

% Steps in the pipeline: 

% 1. Preparing the files
% 2. Segmentation
% 3. Morphometrics analysis
% 4. Foci detection and Cell Cycle analysis

tic
%% 1. Preparing tif files using MIJ (Fiji-Matlab interface)

% Select folder with ND2 files

disp('Select Folder with ND2 files');
dirname = uigetdir();
dirname = fixDir(dirname);


%The macro converts ND2 files to tif for two separate channels. It performs 
%gaussian blur on the green channel

macroFile1='ConvertND2toTif.txt';

% Call MIJ for preprocessing:
filepathMacro = getMacroPath(); %Macro Files path

%variables beginning and end of the movie t end and t benginning

runMacro([filepathMacro,macroFile1],dirname); %calls MIJ to run the Fiji macro without arguments

%% Converting to SuperSegger naming format
%
% The previous script generates files in the following format:
% 01xy - C=0_t0001.tif
% 01xy - C=0_t0002.tif
% ...
% and converts to SuperSegger format:time, xy-positions, fluorescence
% channel, where c1 is bright field or phase contrast  and c2,c3 etc are
% different fluorescent channels
%
% t00001xy001c1.tif
% t00002xy001c1.tif
% etc ... without basename

dirname =  ([dirname,'Analysis',filesep]); %set analysis folder
cd(dirname);
%Set conversion parameters based on your file names:
basename = '';
timeFilterBefore ='_t';
timeFilterAfter = '.t' ;
xyFilterBefore='';
xyFilterAfter='xy';
channelNames = {'C=0','C=1'};

convertImageNames(dirname, basename, timeFilterBefore, ...
    timeFilterAfter, xyFilterBefore,xyFilterAfter, channelNames )

%% Setting the segmentations constants for your bacteria and micrscope resolution
%
% Using correct resolution ensures correct pixel size and segmentation
% constants if you do not know which constants to use you can run
% tryDifferentConstants(dirname) with a phase image to choose. 60X
% indicates 100nm/pix and 100X indicates 60nm/pix

% for E. coli we mainly use : '60XEc' : loadConstants 60X Ecoli - 100
% nm/pix '100XEc': loadConstants 100X Ecoli  - 60 nm/pix

% To see the possible constants type :
%[~, list] = getConstantsList;
% list'

res = '60XBsS750fil2'; %These constants have been optimized for B.subtilis filaments

%% Paralell Processing Mode

% to run code in parallel mode must have the parallel processing toolbox,
% for convenience default is false (non-parallel)

parallel_flag = true;

%% Loading SuperSegger Constants

CONST = loadConstants(res,parallel_flag) ;

%% Modifying specific SuperSegger constant options

% after you load the constants you can modify them according to your needs
% for more options, looks at the loadConstants file.
% The different constant sections are defined below:

%% Channel registration

% Alignment constants for multiple channel files: (Used bead images to
% perform this registration)

% Constants Calarco Microscope:
CONST.imAlign.GFP     = [ 0.0000    0.0000    0.0000    0.0000];
CONST.imAlign.mCherry = [0.04351    -0.0000    0.7200    0.5700]; 
% *last two values are X and Y shifts.

% CONST.imAlign.DAPI    = [-0.0354   -0.0000    1.5500   -0.3900]; these
% are super segger values
% 
% HOW TO OBTAIN REGISTRATION VALUES:
% % im_base = imread('path/of/brightfield/in/non/shifted/channel.tif')
% % im_shifted = imread('path/of/brightfield/in/shifted/channel.tif')
% % [out,~,~] = intAlignIm(im_shifted, im_base, 100 ); %The values in 'out' are the four values needed for the shifted channel.
% 
% Channel order for registration:
CONST.imAlign.out = {CONST.imAlign.GFP,...     % c1 channel (this case is Phase or Bright field)
                    CONST.imAlign.GFP,...      % c2 channel name
                    CONST.imAlign.mCherry,...  % c3 channel name
                    CONST.imAlign.GFP};        % c4 channel name
                
%% Foci detection settings (For use only in SuperSegger)

CONST.trackLoci.numSpots = [5 0]; % Max number of foci to fit in each fluorescence channel (default = [0 0])
CONST.trackLoci.fluorFlag = false ;    % compute integrated fluorescence (default = true)
CONST.trackOpti.NEIGHBOR_FLAG = false; % calculate number of neighbors (default = false)
CONST.imAlign.AlignChannel = 1; % change this if you want the images to be aligned to fluorescence channel
CONST.view.fluorColor = {'b','g','r'}; %Set the color for plotting different channels (in order)

%% Segmentation Filtering options in SuperSegger (modify if you see problems in detection)

%The first two options (PEBBLE_CONST, INTENSITY_DIF) are for the remove_debris variable. 
CONST.superSeggerOpti.PEBBLE_CONST = 1.3; %Default 1.5 for 60XEcM9
CONST.superSeggerOpti.INTENSITY_DIF = 0.3; %Default 0.15 for 60XEcM9
CONST.superSeggerOpti.remove_microcolonies =false; %Default is 1. It prevents deleting clusters of cells.
CONST.superSeggerOpti.remove_debris = 1; %Turn off it is deleting cells
CONST.superSeggerOpti.MAX_WIDTH = 1e15; %Set this high to prevent filaments to be split 
CONST.seg.OPTI_FLAG = false; %To avoid segmenting cells by shape

% this helps cleanup death cells, fragments etc..
CONST.trackOpti.REMOVE_STRAY = true;

% Other constants can be modified using the trainingGui and the sliders,
% specially threshold1 and threshold2 and blur options can improve
% segmentation.
% NOTE: When using trainingGui, Press modify constants twice and load the images to start, 

%% Time step and Pixel size settings:

CONST.getLocusTracks.PixelSize =0.1; %Pixel size in Andor CCD 6.45 um divided by 60X
CONST.getLocusTracks.TimeStep = 5; %time between frames, used for analysis

%% Options to remove CellAsic pillars:

% CONST.superSeggerOpti.remove_pillars.flag = true;
% CONST.superSeggerOpti.remove_pillars.radius = 2;
% CONST.superSeggerOpti.remove_pillars.cut = 0.05;
% CONST.superSeggerOpti.remove_pillars.Area_Cut = 700;
% CONST.superSeggerOpti.remove_pillars.debug = false;

%% Skip Frames for Segmentation

% For fast time-lapse or slow growth you can skip phase image frames 
% during segmentation to increase processing speed. Fluorescence images 
% will not be skipped. 

skip = 1;  % segment every frame
%skip = 5; % segment every fifth phase image

%% Clean previous segmented data

% If set to true, will begin processing on aligned images; if false, will 
% try to restart processing at last successful function (default = false)

cleanflag = 1;

%% 2. Running segmentation in SuperSegger

% Analysis steps in SuperSegger

% 1 : Alignment
% 2 : Segmentation
% 3 : Stripping
% 4 : Linking
% 5 : Cell marker
% 6 : Fluor
% 7 : Foci
% 8 : cellA structrues
% 9 : clist
% 10 : cell files

% in MorphoSegger, only Alignment, Segmentation and Stripping are used and
% indicated in the variable startEnd.

startEnd = [1 3];
BatchSuperSeggerOpti( dirname, skip, cleanflag, CONST,startEnd);

%clearvars -except dirname filepathMacro

%% Clean up files (raw_im *.tif,*.mat)

%Delete original and raw_im folders:
disp('Deleting files...')
rmdir( [dirname,filesep,'original',filesep],'s' ); 
rmdir( [dirname,filesep,'raw_im',filesep],'s' );


%% Converting images to stack 

disp('Converting Images to Stack');
macroFile2='tifImagesToStack.txt';
runMacro([filepathMacro,macroFile2],dirname);

%% 3. Run Morphometrics for Pill Mesh calculation 

% Morphometrics is a robust pipeline that interpolates cell contours at
% subpixel resolution using PSICIC. It uses Oufti routines for segmentation
% and ObjectJ routines for division detection. Here we will use
% Morphometrics to quickly re-segment the binary masks, contour
% fitting, mesh calculation together with basic lineage tracking.

% paramName ='Morphometrics_prefs_mask_CL'; %Select parameter file 
% params = loadParams( paramName );
% 
% %List of most frequently changed parameters, modify here for different
% %types of images
% params.v_imtype = 2;        % 1 = Phase; 2 = Fluorescence (internal); 3 = Fluorescence (peripheral)      
% params.v_method = 3;        % 1 = Gradient Segmentation; 2 = Laplacian Segmentation; 
%                             % 3 = Adaptive Threshold Segmentation; 4 = Canny Segmentation
% params.v_simplethres=1;     % Simple threshold 
% params.f_areamin = 10;     % Min region size
% params.f_areamax = 200000;  % Max region size
% params.v_prox = 0;          % Cells are in proximity
% params.v_exclude=0;         % Exclude edge objects
% params.f_hmin_split=2;      % Cut distance (pxls)
% params.v_save = 1;          % Save output
% params.v_mt_mesh=1;         % Pill mesh
% params.v_falsepos = 1;      % Reject false positives
% params.f_int_rej = 3;       % False positive rejection parameter
% % Tracking parameters:
% params.f_pert_same = 0.55;  % Fractional overlap
% params.f_frame_diff = 4;    % Frame overlap
% %workers = 6;                % Number of workers for parallel job
% 
% disp('Running Morphometrics in parallel')
% run_parallel(dirname,params);

%% 4. Foci calculation - Diego's pipeline

%make cell cycle dir
% loop through 
% parfor

% file_filter = '*.tif';
% dirname_xy = fixDir(dirname_xy);
% 
% % Reset n values in case directories have already been made.
% contents = dir([dirname_xy,'fluor*']);
% num_dir_tmp = numel(contents);
% nc = 1;
% num_c = 1;
% 
% % reset values for nc
% for i = 1:num_dir_tmp
%     if (contents(i).isdir) && (numel(contents(i).name) > numel('fluor'))
%         num_c = num_c+1;
%         nc = [nc, str2double(contents(i).name(numel('fluor')+1:end))+1];
%     end
% end


%% Shutting down parallel pool

disp('Closing parallel pool...')
poolobj = gcp('nocreate');
delete(poolobj);

%% Saving analysis parameters





%% THE END

t1=toc;
disp(['Finished in ' num2str(round(10*t1/60)/10) ' minutes.']);
cd(dirname);

%% Unused code (might be useful for some other things


%disp ('Do you want to cleanup folders?, Y/N [Y]:')
% try
%     answer=input('','s');
%     if lower(answer) ~='y'
%         disp ('Exiting BatchSuperSegger. Reset clean flag and rerun');
%         return
%     end
% catch
%     % can not use input  - in eval mode
% end

%copyfile( [dirname,filesep,'xy1',filesep,'seg',filesep,'*.tif'], [dirname,filesep,'morphometrics'] )


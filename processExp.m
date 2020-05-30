%%%%%%%%%%%%%%%%%%%%% - Pipeline for Imaging Analysis - %%%%%%%%%%%%%%%%%%%
%                                                                         %
%     @@@@@@@@@@@@@@@@@@          @   @        @@@@@@@@@@@@@@@@@@@@       %
%         @@@@@@@@@@@@@@@@@       @@@@@       @@@@@@@@@@@@@@@@@           %
%             @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@              %
%               @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                %
%                 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                 %
%                          @@@@@@@@@@@@@@@@@@@                            %
%                               @@@@@@@@@                                 %
%                                  @@@@                                   %
%                                   @                                     %
%                                                                         %
% Process ND2 files, split channels, converts to tif, performs segmentation
% with SuperSegger and analysis with Morphometrics v2. . The current
% script was modified from SuperSegger. Read full license

% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

%%  MODIFY PIPELINE ACCORDING TO EXPERIMENT SETTINGS AND SAVE IT IN EXPERIMENT FOLDER

function processExp(preprocess,naming,supersegger,cleanup,imag2stack)

% Steps of the pipeline: 
%Boolean variables to decide which steps of the pipeline to run:

% 1. Preparing files           - preprocess
% 2. Convert file names        - naming
% 3. SuperSegger Segmentation  - supersegger
% 4. Clean up  files           - cleanup
% 5. Convert images to stack   - imag2stack
% 6. Morphometrics             - morpho
% 7. Foci analysis             - fociCalc

tic
%% 1. Preparing tif files using MIJ (Fiji-Matlab interface)

if preprocess  
    
    % Select folder with ND2 files
    disp('Select Folder with ND2 files:');
    dirname = uigetdir();
    dirname = fixDir(dirname);

    %The macro converts ND2 files to tif for two separate channels. It performs 
    %gaussian blur on the green channel if desired (check the macro code)

    macroFile1='ConvertND2toTif.txt';

    % Call MIJ for preprocessing:
    filepathMacro = getMacroPath(); %Macro Files path

    % Select frames to analyze
    t_start=1;
    t_end=89;

    %calls MIJ to run the Fiji macro with arguments
    args=strcat(dirname,';',num2str(t_start),';',num2str(t_end)); %group arguments
    runMacro([filepathMacro,macroFile1],args); 

end
%% 2. Converting to SuperSegger naming format
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

if naming  
    
    if ~exist( 'dirname','var')
        dirname= pwd;
        dirname = fixDir(dirname);
        dirname =  ([dirname,'Analysis',filesep]);
    end        
    
    dirname =  ([dirname,'Analysis',filesep]);
    cd(dirname); % **This is important to keep**

    %Set conversion parameters based on your file names:
    basename = '';
    timeFilterBefore ='_t';
    timeFilterAfter = '.t' ;
    xyFilterBefore='';
    xyFilterAfter='xy';
    channelNames = {'C=0','C=1'};

    convertImageNames(dirname, basename, timeFilterBefore, ...
        timeFilterAfter, xyFilterBefore,xyFilterAfter, channelNames )

end
%% 3. Setting segmentation constants and run SuperSegger
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

if supersegger    
    if ~exist( 'dirname','var')
        dirname= pwd;
        dirname = fixDir(dirname);
        dirname =  ([dirname,'Analysis',filesep]);
    end        

%% Paralell Processing Mode

% to run code in parallel mode must have the parallel processing toolbox,
% for convenience default is false (non-parallel)

    parallel_flag = true;

%% Loading SuperSegger Constants

    res = '60XBsS750fil2'; %These constants have been optimized for B.subtilis filaments
    CONST = loadConstants(res,parallel_flag) ;

%% Modifying specific SuperSegger constant options

% after you load the constants you can modify them according to your needs
% for more options, looks at the loadConstants file.
% The different constant sections are defined below:

%% Channel registration

% Alignment constants for multiple channel files: (Used bead images to
% perform this registration)

    % Constants Calarco Microscope:
    CONST.imAlign.Phase     = [ 0.0000    0.0000    0.0000    0.0000];
    CONST.imAlign.GFP       = [ 0.0000    0.0000    0.0000    0.0000];
    CONST.imAlign.mCherry   = [0.04351   -0.0000    0.7200    0.5700]; 
    CONST.imAlign.DAPI      = [0.0000     0.0000    0.0000    0.0000]; 

% *last two values are X and Y shifts.
 
% HOW TO OBTAIN REGISTRATION VALUES:
% % im_base = imread('path/of/brightfield/in/non/shifted/channel.tif')
% % im_shifted = imread('path/of/brightfield/in/shifted/channel.tif')
% % [out,~,~] = intAlignIm(im_shifted, im_base, 100 ); %The values in 'out' are the four values needed for the shifted channel.
% 
% Channel order for registration:
    CONST.imAlign.out = {CONST.imAlign.Phase,...     % c1 channel (this case is Phase or Bright field)
                        CONST.imAlign.GFP,...      % c2 channel name
                        CONST.imAlign.mCherry,...  % c3 channel name
                        CONST.imAlign.DAPI};        % c4 channel name

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
    CONST.superSeggerOpti.MAGIC_RADIUS = 7; % radius of contrast enhancement. deletes areas between cells 
    CONST.seg.OPTI_FLAG = false; %To avoid segmenting cells by shape
    CONST.regionOpti.MIN_LENGTH = 0; % min length of cells
    CONST.trackOpti.MIN_AREA=80; %filter small particles


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

%% Running segmentation in SuperSegger

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

    startEnd = [1 10];
    BatchSuperSeggerOpti( dirname, skip, cleanflag, CONST,startEnd);


end
%% 4. Clean up files (raw_im *.tif, original, *.tif)

if cleanup
     if ~exist( 'dirname','var')
        dirname= pwd;
        dirname = fixDir(dirname);
        dirname =  ([dirname,'Analysis',filesep]);
    end     
    
    %Delete original and raw_im folders:
    disp('Deleting files...')
    delete( [dirname,filesep,'original',filesep,'*.tif' ]);
    delete( [dirname,filesep,'raw_im',filesep,'*.tif' ]);

end
%% 5. Converting images to stack 

if imag2stack   
     if ~exist( 'dirname','var')
        dirname= pwd;
        dirname = fixDir(dirname);
        dirname =  ([dirname,'Analysis',filesep]);
     end     
    
    disp('Converting Images to Stack...');
    macroFile2='tifImagesToStack.txt';
    runMacro([filepathMacro,macroFile2],dirname);

end

%% THE END

t1=toc;
disp(['Finished in ' num2str(round(10*t1/60)/10) ' minutes.']);

load('handel') %alarm that the code is finished
sound(y,Fs)
end



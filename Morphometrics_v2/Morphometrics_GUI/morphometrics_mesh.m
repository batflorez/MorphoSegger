%Tristan Ursell
%March 2012
%
%This function takes a set points defined by X and Y, that are presumed to
%form a closed contour and returns a uniquely-determined mesh and branching
%neutral axis diagram, appended to the input data structure:
%
% frame.object
%
%In particular, this script will add to this structure two sub-fields:
%
% frame.object.mesh
% frame.object.branch
% frame.object.width
%
%The 'mesh' sub-field has the sub-fields:
%
% X_pts = the x points of the current mesh polygon
% Y_pts = the y points of the current mesh polygon
% area = the area of the current mesh polygon
% perim = the perimeter length of the current mesh polygon
% centX = the x coordinate of the current mesh polygon centroid
% centY = the y coordinate of the current mesh polygon centroid
% neighbors = indices of the neighboring mesh polygons
%
%The 'branch' sub-field has the sub-fields:
%
% Xpos = the x coordinates of the current neutral axis
% Ypos = the y coordinates of the current neutral axis
% degree = the degree of connection of each point on the current neutral
% axis. A degree = 1 indicates a branch end-point that does not connect to
% another branch. A degree = 3 or higher indicates a branch end-point that
% meets that number of other branch end-points. A degree = 2 indicates the
% points is located within a branch.
% neighbors_left = the indices of all other branches that connect at one
% end of the current branch.
% neighbors_right = the indices of all other branches that connect at the
% other end of the current branch.
%


function morphometrics_mesh(handles,file1,path1)
global home_path
cd(home_path)

if isempty(file1)
    %load mat data file
    [file1,path1]=uigetfile({'*.mat','mat-data file'},'Select Contour Data File');
    if file1==0
        return
    end
end

%set path
home_path=path1;
cd(home_path)

set(handles.text_process,'String','loading file ...')
drawnow()
load(fullfile(path1,file1));
set(handles.text_process,'String','processing ...')
drawnow()

if ~exist('frame','var')
    disp('This mat file does not contain the "frame" structure required for meshing.')
    return
end

Nframe=length(frame);

%get current contour
for k=1:Nframe
    for c=1:length(frame(k).object)
        
        %check for stop
        global stopV
        if stopV
            stopV=0;
            set(handles.text_process,'String',['Meshing stopped.'])
            return
        end
        
        clearvars -except frame outname v_* f_* cell path1 file1 k c h0 handles eventdata hObject Nframe Ncell
        
        if isfield(frame(k).object(c),'Xcont')
            %get contour points
            X=frame(k).object(c).Xcont';
            Y=frame(k).object(c).Ycont';
            
            %make sure contours are long enough
            if length(X)>4
                %find meshes and branches / possible plotting
                %if get(handles.checkbox_meshplot,'Value')
                %    [frame(k).object(c).mesh,frame(k).object(c).branch]=contour2mesh(X,Y,handles,'plot');
                %else
                [frame(k).object(c).mesh,frame(k).object(c).branch]=contour2mesh(X,Y,handles);
                %end
            end
            set(handles.text_process,'String',['meshing: frame ' num2str(k) ' of ' num2str(Nframe) ', cell ' num2str(frame(k).object(c).cellID)])
            drawnow
        end
    end
end

fname=fullfile(path1,file1);
[~,bname,~]=fileparts(fname);
fname=fullfile(path1,[bname '_full_MESH.mat']);
disp(['Mesh and branch data added to:  ' [bname '_full_MESH.mat']])
set(handles.text_process,'String','Finished meshing.')
clearvars -except frame* cell* outname fname v_* f_* Nframe Ncell
save(fname)







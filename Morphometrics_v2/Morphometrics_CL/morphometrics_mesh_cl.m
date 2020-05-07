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
% Saves cells_ID_table - Andres Florez


function morphometrics_mesh_cl(contourname)

load(contourname);
disp('processing ...')

if ~exist('frame','var')
    disp('This mat file does not contain the "frame" structure required for meshing.')
    return
end

Nframe=length(frame);

%get current contour
for k=1:Nframe
    for c=1:length(frame(k).object)
        clearvars -except frame cells_ID_table outname v_* f_* cell outname k c h0 eventdata hObject Nframe Ncell
        
        if isfield(frame(k).object(c),'Xcont')
            %get contour points
            X=frame(k).object(c).Xcont';
            Y=frame(k).object(c).Ycont';
            
            %make sure contours are long enough
            if length(X)>4
                %find meshes and branches / possible plotting
                %[frame(k).object(c).mesh,frame(k).object(c).branch]=contour2mesh(X,Y,'plot');
                [frame(k).object(c).mesh,frame(k).object(c).branch]=contour2mesh(X,Y);
            end
            disp(['meshing: frame ' num2str(k) ' of ' num2str(Nframe) ', cell ' num2str(frame(k).object(c).cellID)])
        end
    end
end

[path1,bname,~]=fileparts(outname);
fname=fullfile(path1,[bname '_full_MESH.mat']);
disp(['Mesh and branch data added to:  ' [bname '_full_MESH.mat']])
clearvars -except frame* cells_ID_table outname fname v_* f_* Nframe Ncell
save(fname)







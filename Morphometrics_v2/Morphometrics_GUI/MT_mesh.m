function MT_mesh(handles,file1,path1)
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
for k=1:length(frame)
    for c=1:length(frame(k).object)
        clearvars -except frame outname v_* f_* cell path1 file1 k c h0 handles eventdata hObject Nframe Ncell
        
        if isfield(frame(k).object(c),'Xcont')
            %get contour points
            X=frame(k).object(c).Xcont(:);
            Y=frame(k).object(c).Ycont(:);
            kappa=frame(k).object(c).kappa_smooth;
            
            %find meshes
            %mesh1 = contour2pill_mesh([X,Y],0.5,0.1,20);
            
            mesh1 = contour2pill_mesh([X,Y],1,0.01,20);
            if mesh1==-1
                mesh1 = contour2pill_mesh([X,Y],1,0.01,50);
            end
            
            frame(k).object(c).pill_mesh=mesh1;
            
            if and(~isempty(mesh1),mesh1>0)
                [ steplength,...
                    frame(k).object(c).centerline,...
                    frame(k).object(c).length,...
                    frame(k).object(c).width] =...
                    findlength(mesh1);
                
                %{
                %test plot
                if c==2
                    plot(X,Y,'k')
                    hold on
                    plot(frame(k).object(c).centerline(:,1),frame(k).object(c).centerline(:,2),'b')
                    for i=1:size(mesh1,1)
                        plot([mesh1(i,1),mesh1(i,3)],[mesh1(i,2),mesh1(i,4)],'r')
                    end
                    hold off
                    axis equal tight
                    xlabel('X')
                    ylabel('Y')
                    drawnow
                end
                %}
                
                %find various curvatures
                frame(k).object(c).cent_kappa=...
                    contour2curvature(frame(k).object(c).centerline(:,1),frame(k).object(c).centerline(:,2));
                
                %calculate the curvatures of the countour along the
                %centerline normal vector, i.e. the association between
                %centerline curvature and contour curvature.
                
                %contour intersection point 1
                x1=[frame(k).object(c).pill_mesh(:,1)];
                y1=[frame(k).object(c).pill_mesh(:,2)];
                
                %contour intersection point 2
                x2=[frame(k).object(c).pill_mesh(:,3)];
                y2=[frame(k).object(c).pill_mesh(:,4)];
                
                %find two cloest contour points to intersection
                cc1=zeros(size(mesh1,1),1);
                cc2=zeros(size(mesh1,1),1);
                for p=2:sum(size(mesh1,1))-1
                    %first point
                    D1=sqrt((X-x1(p)).^2 + (Y-y1(p)).^2);
                    [d1,i1]=sort(D1);
                    
                    d1_sum=d1(1)+d1(2);
                    cc1(p)=(kappa(i1(1))*d1(2)+kappa(i1(2))*d1(1))/d1_sum;
                    
                    %second point
                    D2=sqrt((X-x2(p)).^2 + (Y-y2(p)).^2);
                    [d2,i2]=sort(D2);
                    
                    d2_sum=d2(1)+d2(2);
                    cc2(p)=(kappa(i2(1))*d2(2)+kappa(i2(2))*d2(1))/d2_sum;
                    
                    %check
                    if or(abs(i1(1)-i1(2))~=1,abs(i2(1)-i2(2))~=1)
                        disp('warning:  possible problem with point association.')
                        disp(['frame: ' num2str(k) ', object: ' num2str(c) ', center line point: ' num2str(p)])
                    end
                end
                
                %record correlated curvature values
                frame(k).object(c).side1_kappa=cc1;
                frame(k).object(c).side2_kappa=cc2;
                
                %calculate the principle curvatures of the projected 2D surface
                %Gaussian curvature (in inverse pixels) is kappa_p1*kappa_p2
                %Mean curvature (in inverse pixels) is 1/2*(kappa_p1+kappa_p2)
                
                maxL=frame(k).object(c).length;
                [frame(k).object(c).kappa_p1,frame(k).object(c).kappa_p2]=...
                    contour2curve_tensor(X,Y,frame(k).object(c).centerline(:,1),frame(k).object(c).centerline(:,2),maxL);
                
                set(handles.text_process,'String',['meshing: frame ' num2str(k) ' of ' num2str(length(frame)) ', cell ' num2str(frame(k).object(c).cellID)])
                disp(['meshing: frame ' num2str(k) ' of ' num2str(length(frame)) ', cell ' num2str(frame(k).object(c).cellID)])
            else
                disp(['Could not create pill mesh for frame ' num2str(k) ', object ' num2str(c) '.'])
                set(handles.text_process,'String',['Could not create pill mesh for frame ' num2str(k) ', object ' num2str(c) '.'])
                
                frame(k).object(c).pill_mesh=[];
                frame(k).object(c).centerline=[];
                frame(k).object(c).length=[];
                frame(k).object(c).width=[];
                frame(k).object(c).cent_kappa=[];
                frame(k).object(c).side1_kappa=[];
                frame(k).object(c).side2_kappa=[];
                frame(k).object(c).kappa_p1=[];
                frame(k).object(c).kappa_p2=[];
            end
        end
        drawnow
    end
end

fname=fullfile(path1,file1);
[~,bname,~]=fileparts(fname);
fname=fullfile(path1,[bname '_pill_MESH.mat']);
disp(['Mesh and length data added to:  ' [bname '_pill_MESH.mat']])
set(handles.text_process,'String','Finished meshing.')
clearvars -except frame* cell* outname fname v_* f_* Nframe Ncell
save(fname)

function [steplength,centerline,length,width] = findlength(mesh)
%calculate centerline
centerline(:,1)=(mesh(:,1)+mesh(:,3))/2;
centerline(:,2)=(mesh(:,2)+mesh(:,4))/2;

%calculate width
width=sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);

%calculate length
x1=(mesh(2:end,1)+mesh(2:end,3))/2;
x2=(mesh(1:end-1,1)+mesh(1:end-1,3))/2;
y1=(mesh(2:end,2)+mesh(2:end,4))/2;
y2=(mesh(1:end-1,2)+mesh(1:end-1,4))/2;

steplength = sqrt((x2-x1).^2+(y2-y1).^2);
length = sum(steplength);




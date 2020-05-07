%Tristan Ursell
%Invert 'frame' and 'cellID' data structures to generate cellID table.
%May 2015

function invert_frame(hObject, eventdata, handles,varargin)

%get file path
if isempty(varargin)
    [filename, pathname] = uigetfile({'*.mat'},'Select mat-file.');
    path_in=fullfile(pathname,filename);
else
    path_in=varargin{1};
end

if filename==0
    return
end

%load mat file
load(path_in)

%check data types
if v_indpt==0
    chk1=1;
else
    chk1=0;
    warning('The cells in this file may have been processed as separate objects.')
end

if isfield(frame(1).object,'cellID')
    chk2=1;
else
    chk2=0;
    error('Could not find the necessary data to construct the cell table.')
end

%construct table
for i=1:Ncell
    set(handles.text_process,'String',['Processing cell: ' num2str(i)])
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

%save table
save(path_in,'cells','-append')

set(handles.text_process,'String',['Cell tabled created.'])

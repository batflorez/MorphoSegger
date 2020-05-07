function [cell,trackframe]=manual_tracker(basename,fname,Nim,handles)
%Tristan Ursell
%May 2014
%
%
%cell(k).bw_label = all of the bw_labels assigned to cell 'k' in the frames given by cell(k).frames
%

%declare global variables
global v_imtype

%initiate 'cellID' field in the 'frame' data structure
%[trackframe().cellID]=0;
for j=1:Nim
    trackframe(j).cellID=0;
end

%get image size
info1=imfinfo(fname);
sx=info1.Width;
sy=info1.Height;

%constructing the lineage (which cells are which in each frame)
f1=figure('units','normalized','outerposition',[0 0 1 1]);
break_opt=0;
Nmax=0;
%loop over actual cell identities
while break_opt==0
    %set index variables
    Nmax=Nmax+1;
    j=0;
    while j<Nim
        %index
        j=j+1;
        if j<1
            j=1;
        end
        
        if j==1
            %load raw and segmented images
            if v_imtype==1
                raw_j=mat2gray(-double(imread(fname,j)));
            else
                raw_j=mat2gray(double(imread(fname,j)));
            end
            Lfinal_j=imread([basename '_Gparent.tif'],j);
            
            %identify assigned lineages
            already_assigned=find(trackframe(j).cellID~=0);
            
            %construct assigned lineage mask
            already_mask=zeros(size(raw_j));
            for m=1:length(already_assigned)
                already_mask = already_mask + (Lfinal_j==already_assigned(m));
            end
            
            %show first image for inital cell selection
            hold off
            view1(:,:,1)=raw_j.*(Lfinal_j==0);
            view1(:,:,2)=raw_j.*~already_mask;
            view1(:,:,3)=raw_j.*((Lfinal_j==0) + already_mask);
            
            clf
            imagesc(view1)
            hold on
            axis equal tight
            title(['Frame: ' num2str(j) ', Select cell: ' num2str(Nmax)]) %assigns frame nr. num2str. converts number into text.
            box on
            xlabel('X')
            ylabel('Y')
        else
            %load raw and segmented images
            if v_imtype==1
                raw_jm1=mat2gray(-double(imread(fname,j-1)));
                raw_j=mat2gray(-double(imread(fname,j)));
            else
                raw_jm1=mat2gray(double(imread(fname,j-1)));
                raw_j=mat2gray(double(imread(fname,j)));
            end
            Lfinal_j=imread([basename '_Gparent.tif'],j);
            Lfinal_jm1=imread([basename '_Gparent.tif'],j-1);
            
            %find the object number of the previousy clicked cell
            curr_obj=find([trackframe(j-1).cellID]==Nmax);
            
            curr_im=zeros(size(raw_j));
            if ~isempty(curr_obj)
                for m=1:length(curr_obj)
                    curr_im=curr_im + bwperim(Lfinal_jm1==curr_obj(m));
                end
            end
            
            %identify assigned lineages
            already_assigned=find(trackframe(j).cellID~=0);
            
            %construct assigned lineage mask
            already_mask=zeros(size(raw_j));
            for m=1:length(already_assigned)
                already_mask = already_mask + (Lfinal_j==already_assigned(m));
            end
            
            %show first image for inital cell selection
            view1(:,:,1)=raw_j.*~curr_im;
            view1(:,:,2)=raw_j.*(Lfinal_jm1==0);
            view1(:,:,3)=raw_j.*(Lfinal_jm1==0) + curr_im;
            
            view2(:,:,1)=raw_j.*(Lfinal_j==0);
            view2(:,:,2)=raw_j.*~already_mask;
            view2(:,:,3)=raw_j.*((Lfinal_j==0) + already_mask);
            
            clf
            imagesc(0.33*view1+0.67*view2)
            axis equal tight
            title(['Frame: ' num2str(j) ', Select cell: ' num2str(Nmax) '  (click outside image for help)']) %assigns frame nr. num2str. converts number into text.
            box on
            xlabel('X')
            ylabel('Y')
        end
        
        %directions string, does not include 'combine' command
        directions_string={'Object color meanings:','objects in current frame = green','objects in previous frame = red',...
            'previously selected objects = blue outline','objects from previous lineage = blue fill','unsegmented regions = grayscale',...
            ' ','Key strokes:','left arrow = go back one frame','right arrow = go forward one frame',...
            'space bar = finish current lineage, start new one','delete = stop all lineage tracking, return to Morphometrics'};
        
        %{
        %directions string, does include 'combine' command
        directions_string={'Colors:','current frame = green','previous frame = red',...
            'previous selection = blue outline','previous lineage = blue','unsegmented regions = grayscale',...
            ' ','Key strokes:','left arrow = go back one frame','right arrow = go forward one frame',...
            'c = combine selected regions','b = branch current lineage','space bar = halt current lineage',...
            'delete = stop all lineage tracking'};
        %}
        
        %{
        annotation(gcf,'rectangle',[0.75 0.53 0.125 0.22],'Linewidth',2);
        
        annotation(gcf,'textbox',...
            [0.75 0.5 0.15 0.25],...
            'String',directions_string,...
            'FitBoxToText','off','Linestyle','none');
        %}        
                
        %select point in current region of interest
        [curr_cell_x,curr_cell_y,button1]=ginput(1);
        
        %if user clicks on non-region, give option list
        if sum(button1==[29,28,32,99,8,98])>0
            %right arrow = goes forward a frame and begins lineage there
            %back arrow = goes back one frame, and allows for re-selection (may be used in tandem)
            %space bar = halts the current cell lineage, and moves onto the next cell lineage
            %delete/backspace = halt lineage algorithm entirely
            %c = combine more than one region into a single cell in the current frame
            
            %right arrow
            if button1==29
                %title(['Frame: ' num2str(j) ' (skipped) , Select cell: ' num2str(Nmax)]) %assigns frame nr. num2str. converts number into text.
                %pause(0.2)
            end
            
            %left arrow
            if and(j>1,button1==28) 
                temp1=find([trackframe(j-1).cellID]==Nmax);
                trackframe(j-1).cellID(temp1)=0;
                j=j-2;
            elseif and(j==1,button1==28)
                temp1=find([trackframe(j).cellID]==Nmax);
                trackframe(j).cellID(temp1)=0;
                j=j-1;
            end
            
            if button1==32 %space bar
                j=Nim;
            end
            
            if button1==99 %c (combine)
                inval1=1;
                while inval1
                    title(['Frame: ' num2str(j) ', Select cell: ' num2str(Nmax) ' (hit enter to finish)']) %assigns frame nr. num2str. converts number into text.
                    %select point in current region of interest
                    [curr_cell_x,curr_cell_y]=getpts;
                    
                    %fix out of bound values
                    curr_cell_x(curr_cell_x<1)=1;
                    curr_cell_y(curr_cell_y<1)=1;
                    curr_cell_x(curr_cell_x>sx)=sx;
                    curr_cell_y(curr_cell_y>sy)=sy;
                    
                    %determine which region was clicked
                    coordx=round(curr_cell_x);
                    coordy=round(curr_cell_y);
                    
                    curr_region=zeros(length(coordx),1);
                    for m=1:length(coordx)
                        curr_region(m)=Lfinal_j(coordy(m),coordx(m));
                    end
                    
                    %check for invalid regions
                    if sum(curr_region==0)>0
                        title('Invalid region(s) selected, try again.')
                        pause(1)
                    else
                        inval1=0;
                    end
                end
                %save the cellID into the frame data
                trackframe(j).cellID(curr_region(curr_region>0))=Nmax;
            end
            
            %delete
            if button1==8 
                break_opt=1;
                break
            end
            
            %handle lineage breaking
            if button1==98
            end
            
        elseif button1==1
            %if external region selected, notify
            if any([curr_cell_x<1,curr_cell_y<1,curr_cell_x>sx,curr_cell_y>sy])
                title('Invalid region selected, try again.')
                uiwait(msgbox(directions_string,'Manual Tracker Directions','modal'))
                j=j-1;
            else
                %determine which region was clicked
                coordx=round(curr_cell_x);
                coordy=round(curr_cell_y);
                curr_region=Lfinal_j(coordy,coordx);
                
                %catch clicks not in a region
                if sum(curr_region==0)>0
                    title('Invalid region selected, try again.')
                    pause(1)
                    j=j-1;
                else
                    %save the cellID into the frame data
                    trackframe(j).cellID(curr_region)=Nmax;
                end
            end
        else
            title('Invalid key stroke.')
            pause(1)
            j=j-1;
        end
        clear curr_region
    end
    %if 'stop all' selected, halt lineage tracking entirely
    if break_opt
        break
    end
end

clear cell
%initialize structure output
for k=1:max([trackframe.cellID])
    N_ents=sum([trackframe.cellID]==k);
    cell(k).frames=zeros(1,N_ents);
    cell(k).bw_label=zeros(1,N_ents);
end

%create cell output with same form as the automatic tracking algorithm
for j=1:Nim
    %find all bwlabels in current frame
    bw_temp=find(trackframe(j).cellID~=0);
    
    if isempty(bw_temp)
        continue
    end
    
    %find all current cellIDs in current frame
    ID_temp=trackframe(j).cellID(bw_temp);
    
    %populate the cell array
    for k=1:length(ID_temp)
        mark1=find(cell(ID_temp(k)).frames==0,1,'first');
        
        cell(ID_temp(k)).frames(mark1)=j;
        cell(ID_temp(k)).bw_label(mark1)=bw_temp(k);
        
        %cell(ID_temp(k)).frames=[cell(ID_temp(k)).frames,j];
        %cell(ID_temp(k)).bw_label=[cell(ID_temp(k)).bw_label,bw_temp(k)];
    end
end

close(f1)

